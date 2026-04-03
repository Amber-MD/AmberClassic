#include "NonBulkGPU.h"
#include <cassert>
#include <thrust/complex.h>
#include <thrust/device_vector.h>
#include <cstdio>

namespace {

  const int WARP_SIZE = 32;
  const int MAX_N_SCATTER_TYPES = 16;

  /***
   * Alternative implementation to calc_f_non_bulk_kernel
   *    with manual GPU cache management
   *
   * Original code taken from
   *    kXray.cu:kXrayGetDerivative1_kernel() [by David S. Cerutti]
   * */
  template<int BLOCK_SIZE, typename FloatType>
  __global__
  void calc_f_non_bulk_w_manual_caching_kernel(int n_atom,
                                               const FloatType* global_frac_xyz,
                                               const FloatType* global_b_factor,
                                               const FloatType* global_occupancy,
                                               int n_hkl,
                                               const int* global_hkl,
                                               const FloatType* global_mss4,
                                               int n_scatter_types,
                                               const FloatType* global_atomic_scatter_factor,
                                               const int* global_scatter_type_index,
                                               thrust::complex<FloatType>* global_f_non_bulk) {
    assert (BLOCK_SIZE == blockDim.x);
    assert (BLOCK_SIZE >= WARP_SIZE);
    assert (BLOCK_SIZE % WARP_SIZE == 0);

    int tgx = threadIdx.x % WARP_SIZE;
    int warp_idx = threadIdx.x / WARP_SIZE;
    int n_wraps = blockDim.x / WARP_SIZE;

    assert(n_wraps > 5); // Needed for parallel reflex data caching (see switch-case below)
    assert(n_scatter_types <= MAX_N_SCATTER_TYPES);

    __shared__ FloatType frac_x[BLOCK_SIZE];
    __shared__ FloatType frac_y[BLOCK_SIZE];
    __shared__ FloatType frac_z[BLOCK_SIZE];
    __shared__ FloatType b_factor[BLOCK_SIZE];
    __shared__ FloatType occupancy[BLOCK_SIZE];
    __shared__ int scatter_type_index[BLOCK_SIZE];

    __shared__ int h[WARP_SIZE];
    __shared__ int k[WARP_SIZE];
    __shared__ int l[WARP_SIZE];
    __shared__ FloatType mss4[WARP_SIZE];
    __shared__ FloatType sf_real[WARP_SIZE];
    __shared__ FloatType sf_imag[WARP_SIZE];

    // indexed by WARP_SIZE * <type> + <hkl_idx>
    __shared__ FloatType atomic_scatter_factor[WARP_SIZE * MAX_N_SCATTER_TYPES];

    for (unsigned hkl_offset = blockIdx.x * WARP_SIZE; hkl_offset < n_hkl; hkl_offset += gridDim.x * WARP_SIZE) {

      { // Cache reflex data
        unsigned i_hkl = hkl_offset + tgx;
        if (i_hkl < n_hkl) {
          switch (warp_idx) {
            case (0):
              mss4[tgx] = global_mss4[i_hkl];
              break;
            case (1):
              h[tgx] = global_hkl[i_hkl * 3];
              break;
            case (2):
              k[tgx] = global_hkl[i_hkl * 3 + 1];
              break;
            case (3):
              l[tgx] = global_hkl[i_hkl * 3 + 2];
              break;
            case (4):
              sf_real[tgx] = 0;
              sf_imag[tgx] = 0;
              break;
            case (5): {
              for (int scatter_type = 0; scatter_type < n_scatter_types; scatter_type++) {
                atomic_scatter_factor[scatter_type * WARP_SIZE + tgx] = \
            global_atomic_scatter_factor[scatter_type * n_hkl + i_hkl];
              }
              break;
            }
            default:
              break;
          }
        }
      }

      __syncthreads();
      // Loop over all atoms
      for (unsigned global_atom_offset = 0; global_atom_offset < n_atom; global_atom_offset += BLOCK_SIZE) {
        // We need to loop over offsets (not indices itself) ALL threads in block execute __syncthreads();

        // Cache atomic data
        {
          unsigned global_i_atom = global_atom_offset + threadIdx.x;
          if (global_i_atom < n_atom) {
            scatter_type_index[threadIdx.x] = global_scatter_type_index[global_i_atom] - 1;
            b_factor[threadIdx.x] = global_b_factor[global_i_atom];
            occupancy[threadIdx.x] = global_occupancy[global_i_atom];
            frac_x[threadIdx.x] = global_frac_xyz[global_i_atom * 3 + 0];
            frac_y[threadIdx.x] = global_frac_xyz[global_i_atom * 3 + 1];
            frac_z[threadIdx.x] = global_frac_xyz[global_i_atom * 3 + 2];
          }
          __syncthreads();
        }

        // Computing f with manual cache is here:

        for (unsigned hkl_cached_idx = warp_idx;
             hkl_cached_idx < min(n_hkl - hkl_offset, WARP_SIZE); hkl_cached_idx += n_wraps) {
          FloatType f_real = 0;
          FloatType f_imag = 0;
          {
            FloatType t_mss4 = mss4[hkl_cached_idx];
            int t_h = h[hkl_cached_idx];
            int t_k = k[hkl_cached_idx];
            int t_l = l[hkl_cached_idx];

            for (unsigned atom_cache_idx = tgx;
                 atom_cache_idx < min(n_atom - global_atom_offset, BLOCK_SIZE); atom_cache_idx += WARP_SIZE) {
              FloatType angle = 2 * M_PI * (
                t_h * frac_x[atom_cache_idx] +
                t_k * frac_y[atom_cache_idx] +
                t_l * frac_z[atom_cache_idx]
              );
              FloatType f = exp(t_mss4 * b_factor[atom_cache_idx]) *
                            atomic_scatter_factor[(scatter_type_index[atom_cache_idx] * WARP_SIZE) + hkl_cached_idx]
                            * occupancy[atom_cache_idx];

              f_real += f * cos(angle);
              f_imag += f * sin(angle);
            }
          }

          const unsigned int warp_mask = 0xffffffff;

          f_real += __shfl_down_sync(warp_mask, f_real, 16);
          f_real += __shfl_down_sync(warp_mask, f_real, 8);
          f_real += __shfl_down_sync(warp_mask, f_real, 4);
          f_real += __shfl_down_sync(warp_mask, f_real, 2);
          f_real += __shfl_down_sync(warp_mask, f_real, 1);

          f_imag += __shfl_down_sync(warp_mask, f_imag, 16);
          f_imag += __shfl_down_sync(warp_mask, f_imag, 8);
          f_imag += __shfl_down_sync(warp_mask, f_imag, 4);
          f_imag += __shfl_down_sync(warp_mask, f_imag, 2);
          f_imag += __shfl_down_sync(warp_mask, f_imag, 1);

          if (tgx == 0) {
            sf_real[hkl_cached_idx] += f_real;
            sf_imag[hkl_cached_idx] += f_imag;
          }

        }
        __syncthreads();
      }

      __syncthreads();
      {
        unsigned i_hkl = hkl_offset + threadIdx.x;
        if (i_hkl < n_hkl && threadIdx.x < WARP_SIZE) {
          FloatType f_real = sf_real[threadIdx.x];
          FloatType f_imag = sf_imag[threadIdx.x];
          global_f_non_bulk[i_hkl] = thrust::complex<FloatType>(f_real, f_imag);
        }
      }
      __syncthreads();
    }
  }

  template<int BLOCK_SIZE, typename FloatType>
  __global__
  void calc_f_non_bulk_kernel(int n_atom,
                              const FloatType* frac_xyz,
                              const FloatType* b_factor,
                              const FloatType* occupancy,
                              int n_hkl,
                              const int* hkl,
                              const FloatType* mss4,
                              int /*n_scatter_types*/,
                              const FloatType* atomic_scatter_factor,
                              const int* scatter_type_index,
                              thrust::complex<FloatType>* f_non_bulk) {
    assert(BLOCK_SIZE == blockDim.x);
    const int tid = threadIdx.x;
    const int i_hkl = blockIdx.x;
    __shared__ thrust::complex<FloatType> term[BLOCK_SIZE];
    term[tid] = {};

    // Basic code to compute f is here; (is it really what is used?)

    if (i_hkl < n_hkl) {
      const FloatType hkl_mss4 = mss4[i_hkl];
      const int h = hkl[i_hkl * 3 + 0];
      const int k = hkl[i_hkl * 3 + 1];
      const int l = hkl[i_hkl * 3 + 2];

      for (int j_atom = tid; j_atom < n_atom; j_atom += BLOCK_SIZE) {
        const FloatType f = std::exp(hkl_mss4 * b_factor[j_atom]) *
                            atomic_scatter_factor[(scatter_type_index[j_atom] - 1) * n_hkl + i_hkl];
        const FloatType angle = 2 * M_PI * (
          frac_xyz[j_atom * 3 + 0] * h +
          frac_xyz[j_atom * 3 + 1] * k +
          frac_xyz[j_atom * 3 + 2] * l
        );
        term[tid] += thrust::complex<FloatType>{f * std::cos(angle), f * std::sin(angle)} * occupancy[j_atom];
      }

      __syncthreads();

      for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1) {
        if (tid < s) {
          term[tid] += term[tid + s];
        }
        __syncthreads();
      }

      if (tid == 0 && i_hkl < n_hkl) {
        f_non_bulk[i_hkl] = term[0];
      }
    }
  }
}

template<xray::NonBulkKernelVersion KERNEL_VERSION, xray::KernelPrecision PRECISION>
xray::NonBulkGPU<KERNEL_VERSION, PRECISION>::NonBulkGPU(
  int n_hkl, const int* hkl, complex_double* f_non_bulk, const double* mSS4, int n_atom,
  const double* b_factor, const double* occupancy, int n_scatter_types, const int* scatter_type_index,
  const double* atomic_scatter_factor) : NonBulk(n_hkl, hkl, f_non_bulk, mSS4, n_atom,
                                                 b_factor, occupancy, n_scatter_types,
                                                 scatter_type_index, atomic_scatter_factor) {
  m_dev_frac_xyz = thrust::device_vector<FloatType>(n_atom * 3);
  m_dev_b_factor = thrust::device_vector<FloatType>(m_b_factor, m_b_factor + n_atom);
  m_dev_occupancy = thrust::device_vector<FloatType>(m_occupancy, m_occupancy + n_atom);
  m_dev_hkl = thrust::device_vector<int>(m_hkl, m_hkl + m_n_hkl * 3);
  m_dev_mSS4 = thrust::device_vector<FloatType>(m_mSS4, m_mSS4 + m_n_hkl);
  m_dev_atomic_scatter_factor = thrust::device_vector<FloatType>(m_atomic_scatter_factor,
                                                              m_atomic_scatter_factor + m_n_hkl * m_n_scatter_types);
  m_dev_scatter_type_index = thrust::device_vector<int>(m_scatter_type_index, m_scatter_type_index + n_atom);
  m_dev_f_non_bulk = thrust::device_vector<thrust::complex<FloatType>>(m_n_hkl);
}

template<xray::NonBulkKernelVersion KERNEL_VERSION, xray::KernelPrecision PRECISION>
void xray::NonBulkGPU<KERNEL_VERSION, PRECISION>::calc_f_non_bulk(int n_atom, const double* frac_xyz) {
  assert(n_atom == m_n_atom);

  thrust::copy(frac_xyz, frac_xyz + n_atom * 3, m_dev_frac_xyz.begin());
  thrust::fill(m_dev_f_non_bulk.begin(), m_dev_f_non_bulk.end(), 0.0);

  switch (KERNEL_VERSION) {
    case (NonBulkKernelVersion::ManualCaching): {

      const int block_size = 256;
      dim3 numBlocks((m_n_hkl + WARP_SIZE - 1) / WARP_SIZE * WARP_SIZE);
      dim3 threadsPerBlock(block_size);

      calc_f_non_bulk_w_manual_caching_kernel<block_size>
      <<<numBlocks, threadsPerBlock>>>(
        n_atom,
        m_dev_frac_xyz.data().get(),
        m_dev_b_factor.data().get(),
        m_dev_occupancy.data().get(),
        m_n_hkl,
        m_dev_hkl.data().get(),
        m_dev_mSS4.data().get(),
        m_n_scatter_types,
        m_dev_atomic_scatter_factor.data().get(),
        m_dev_scatter_type_index.data().get(),
        m_dev_f_non_bulk.data().get()
      );
      break;
    }
    case (NonBulkKernelVersion::StraightForward): {

      dim3 numBlocks(m_n_hkl);
      const int block_size = 32;
      dim3 threadsPerBlock(block_size);

      calc_f_non_bulk_kernel<block_size>
      <<<numBlocks, threadsPerBlock>>>(
        n_atom,
        m_dev_frac_xyz.data().get(),
        m_dev_b_factor.data().get(),
        m_dev_occupancy.data().get(),
        m_n_hkl,
        m_dev_hkl.data().get(),
        m_dev_mSS4.data().get(),
        m_n_scatter_types,
        m_dev_atomic_scatter_factor.data().get(),
        m_dev_scatter_type_index.data().get(),
        m_dev_f_non_bulk.data().get()
      );
      break;
    }
  }
  thrust::copy(m_dev_f_non_bulk.begin(), m_dev_f_non_bulk.end(), m_f_non_bulk);
}

template class xray::NonBulkGPU<xray::NonBulkKernelVersion::ManualCaching, xray::KernelPrecision::Single>;
template class xray::NonBulkGPU<xray::NonBulkKernelVersion::ManualCaching, xray::KernelPrecision::Double>;
template class xray::NonBulkGPU<xray::NonBulkKernelVersion::StraightForward, xray::KernelPrecision::Single>;
template class xray::NonBulkGPU<xray::NonBulkKernelVersion::StraightForward, xray::KernelPrecision::Double>;
