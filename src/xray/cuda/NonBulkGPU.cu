#include "NonBulkGPU.h"
#include <cassert>
#include <thrust/complex.h>
#include <thrust/device_vector.h>
#include <cstdio>

namespace {

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

    // Basic code to compute f is here; (per xray_non_bulk.cu, this
    // StraightForward code seems to be what is used

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

#if 1  /* first test of symmetrization:  */
      // test for spacegroup_number here; or hard-wire, as here

        const int h2 = -hkl[i_hkl * 3 + 0];
        const int k2 =  hkl[i_hkl * 3 + 1];
        const int l2 = -hkl[i_hkl * 3 + 2];

        const FloatType angle2 = 2 * M_PI * (
          frac_xyz[j_atom * 3 + 0] * h2 +
          frac_xyz[j_atom * 3 + 1] * k2 +
          frac_xyz[j_atom * 3 + 2] * l2
        );
      
        // hard-wire nb=3
        if( k2/3 % 2 != 0 ) {  //N.B.: need na,nb,nc here somehow
          term[tid] -= thrust::complex<FloatType>{f * std::cos(angle2), f * std::sin(angle2)} * occupancy[j_atom];
        } else {
          term[tid] += thrust::complex<FloatType>{f * std::cos(angle2), f * std::sin(angle2)} * occupancy[j_atom];
        }
#endif

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
