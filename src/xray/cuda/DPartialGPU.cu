#include <cassert>

#include "DPartialGPU.h"
#include <thrust/device_vector.h>
#include <thrust/complex.h>
#include <thrust/functional.h>
#include <thrust/transform.h>

#undef P212121
#undef P21
#undef P6
#define NA   1
#define NB   1
#define NC   1

namespace {
  template<int blockDimX, typename FloatType>
  __global__
  void calc_d_target_d_frac_kernel(
    int n_atom,
    const FloatType* frac_by_2_pi,
    const FloatType* b_factor,
    const FloatType* occupancy,
    int n_hkl,
    const int* hkl,
    const FloatType* mSS4,
    const FloatType* atomic_scatter_factor,
    const int* scatter_type_index,
    const FloatType* f_calc_phase,
    const FloatType* f_scale,
    const FloatType* abs_f_calc,
    const FloatType* d_target_d_abs_f_calc_by_2_pi,
    FloatType* d_target_d_frac) {

    const int tid = threadIdx.x;
    const int i = blockIdx.x;
    __shared__ FloatType term_x[blockDimX];
    __shared__ FloatType term_y[blockDimX];
    __shared__ FloatType term_z[blockDimX];
    term_x[tid] = 0;
    term_y[tid] = 0;
    term_z[tid] = 0;

    if (i < n_atom) {

      const int offset = (scatter_type_index[i] - 1) * n_hkl;
      const FloatType atom_b_factor = b_factor[i];
      const FloatType atom_occupancy = occupancy[i];

      for (int i_hkl = tid; i_hkl < n_hkl; i_hkl += blockDimX) {
        if (abs_f_calc[i_hkl] < 1e-3) {
          continue;
        }

        FloatType phase = -(
          hkl[i_hkl * 3 + 0] * frac_by_2_pi[i * 3 + 0] +
          hkl[i_hkl * 3 + 1] * frac_by_2_pi[i * 3 + 1] +
          hkl[i_hkl * 3 + 2] * frac_by_2_pi[i * 3 + 2]
        );

        FloatType scatter_factor = atomic_scatter_factor[offset + i_hkl];

        FloatType f = scatter_factor * exp(mSS4[i_hkl] * atom_b_factor) *
                      f_scale[i_hkl] * d_target_d_abs_f_calc_by_2_pi[i_hkl] *
                      atom_occupancy;
        FloatType tmp = f * sin(phase + f_calc_phase[i_hkl]);

        term_x[tid] += hkl[i_hkl * 3 + 0] * tmp;
        term_y[tid] += hkl[i_hkl * 3 + 1] * tmp;
        term_z[tid] += hkl[i_hkl * 3 + 2] * tmp;

#ifdef P21
        // spacegroup 4 code here

        // set #2:   -h,k,-l
        FloatType phase2 = -(
          - hkl[i_hkl * 3 + 0] * frac_by_2_pi[i * 3 + 0]
          + hkl[i_hkl * 3 + 1] * frac_by_2_pi[i * 3 + 1]
          - hkl[i_hkl * 3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp2 = f * sin(phase2 + f_calc_phase[i_hkl]);
        // note that NB is def-ed at top of file
        if( hkl[i_hkl * 3 + 1]/NB % 2 != 0 ){ tmp2 = -tmp2; }
        term_x[tid] -= hkl[i_hkl * 3 + 0] * tmp2;
        term_y[tid] += hkl[i_hkl * 3 + 1] * tmp2;
        term_z[tid] -= hkl[i_hkl * 3 + 2] * tmp2;
#endif

#ifdef P6
        // spacegroup 168 code here

        // set #2:   h+k,-h,l
        FloatType phase2 = -(
           (hkl[i_hkl*3 + 0] + hkl[i_hkl*3 + 1]) * frac_by_2_pi[i * 3 + 0]
          - hkl[i_hkl*3 + 0] * frac_by_2_pi[i * 3 + 1]
          + hkl[i_hkl*3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp2 = f * sin(phase2 + f_calc_phase[i_hkl]);
        term_x[tid] += (hkl[i_hkl*3 + 0] + hkl[i_hkl*3 + 1])  * tmp2;
        term_y[tid] -= hkl[i_hkl*3 + 0] * tmp2;
        term_z[tid] += hkl[i_hkl*3 + 2] * tmp2;

        // set #3:   k,-h-k,l
        FloatType phase3 = -(
           (hkl[i_hkl*3 + 0] + hkl[i_hkl*3 + 1]) * frac_by_2_pi[i * 3 + 0]
          - hkl[i_hkl*3 + 0] * frac_by_2_pi[i * 3 + 1]
          + hkl[i_hkl*3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp3 = f * sin(phase3 + f_calc_phase[i_hkl]);
        term_x[tid] += hkl[i_hkl*3 + 1] * tmp3;
        term_y[tid] -= (hkl[i_hkl*3 + 0] + hkl[i_hkl*3 + 1])  * tmp3;
        term_z[tid] += hkl[i_hkl*3 + 2] * tmp3;

        // set #4:   -h,-k,l
        FloatType phase4 = -(
          - hkl[i_hkl*3 + 0] * frac_by_2_pi[i * 3 + 0]
          - hkl[i_hkl*3 + 1] * frac_by_2_pi[i * 3 + 1]
          + hkl[i_hkl*3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp4 = f * sin(phase4 + f_calc_phase[i_hkl]);
        term_x[tid] -= hkl[i_hkl*3 + 0] * tmp4;
        term_y[tid] -= hkl[i_hkl*3 + 1] * tmp4;
        term_z[tid] += hkl[i_hkl*3 + 2] * tmp4;

        // set #5:   -h-k,h,l
        FloatType phase5 = -(
          -(hkl[i_hkl*3 + 0] + hkl[i_hkl*3 + 1]) * frac_by_2_pi[i * 3 + 0]
          + hkl[i_hkl*3 + 0] * frac_by_2_pi[i * 3 + 1]
          + hkl[i_hkl*3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp5 = f * sin(phase5 + f_calc_phase[i_hkl]);
        term_x[tid] -= (hkl[i_hkl*3 + 0] + hkl[i_hkl*3 + 1])  * tmp5;
        term_y[tid] += hkl[i_hkl*3 + 0] * tmp5;
        term_z[tid] += hkl[i_hkl*3 + 2] * tmp5;

        // set #6:   -k,h+k,l
        FloatType phase6 = -(
          - hkl[i_hkl*3 + 1] * frac_by_2_pi[i * 3 + 0]
          + (hkl[i_hkl*3 + 0] +hkl[i_hkl*3 + 1]) * frac_by_2_pi[i * 3 + 1]
          + hkl[i_hkl*3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp6 = f * sin(phase6 + f_calc_phase[i_hkl]);
        term_x[tid] -= hkl[i_hkl*3 + 1] * tmp6;
        term_y[tid] += (hkl[i_hkl*3 + 0] + hkl[i_hkl*3 + 1])  * tmp6;
        term_z[tid] += hkl[i_hkl*3 + 2] * tmp6;

#endif

#ifdef P212121
        // spacegroup 19 code here

        // set #2: -h,-k,l
        FloatType phase2 = -(
          - hkl[i_hkl * 3 + 0] * frac_by_2_pi[i * 3 + 0]
          - hkl[i_hkl * 3 + 1] * frac_by_2_pi[i * 3 + 1]
          + hkl[i_hkl * 3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp2 = f * sin(phase2 + f_calc_phase[i_hkl]);
        if( (hkl[i_hkl * 3 + 0]/NA + hkl[i_hkl * 3 + 2]/NC) % 2 != 0 ){ tmp2 = -tmp2; }
        term_x[tid] -= hkl[i_hkl * 3 + 0] * tmp2;
        term_y[tid] -= hkl[i_hkl * 3 + 1] * tmp2;
        term_z[tid] += hkl[i_hkl * 3 + 2] * tmp2;

        // set #3: -h,k,-l
        FloatType phase3 = -(
          - hkl[i_hkl * 3 + 0] * frac_by_2_pi[i * 3 + 0]
          + hkl[i_hkl * 3 + 1] * frac_by_2_pi[i * 3 + 1]
          - hkl[i_hkl * 3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp3 = f * sin(phase3 + f_calc_phase[i_hkl]);
        if( (hkl[i_hkl * 3 + 1]/NB + hkl[i_hkl * 3 + 2]/NC) % 2 != 0 ){ tmp3 = -tmp3; }
        term_x[tid] -= hkl[i_hkl * 3 + 0] * tmp3;
        term_y[tid] += hkl[i_hkl * 3 + 1] * tmp3;
        term_z[tid] -= hkl[i_hkl * 3 + 2] * tmp3;

        // set #4:  h,-k,-l
        FloatType phase4 = -(
            hkl[i_hkl * 3 + 0] * frac_by_2_pi[i * 3 + 0]
          - hkl[i_hkl * 3 + 1] * frac_by_2_pi[i * 3 + 1]
          - hkl[i_hkl * 3 + 2] * frac_by_2_pi[i * 3 + 2]
        );
        FloatType tmp4 = f * sin(phase4 + f_calc_phase[i_hkl]);
        if( (hkl[i_hkl * 3 + 0]/NA + hkl[i_hkl * 3 + 1]/NB) % 2 != 0 ){ tmp4 = -tmp4; }
        term_x[tid] += hkl[i_hkl * 3 + 0] * tmp4;
        term_y[tid] -= hkl[i_hkl * 3 + 1] * tmp4;
        term_z[tid] -= hkl[i_hkl * 3 + 2] * tmp4;

#endif
      }

      __syncthreads();

      for (int s = blockDimX / 2; s > 0; s >>= 1) {
        if (tid < s) {
          term_x[tid] += term_x[tid + s];
          term_y[tid] += term_y[tid + s];
          term_z[tid] += term_z[tid + s];
        }
        __syncthreads();
      }

      if (tid == 0) {
        d_target_d_frac[i * 3 + 0] = term_x[0];
        d_target_d_frac[i * 3 + 1] = term_y[0];
        d_target_d_frac[i * 3 + 2] = term_z[0];
      }
    }
  }

  template<typename T>
  void multiply_by(T factor, thrust::device_vector<T> &dev_vec) {
    using namespace thrust::placeholders;
    thrust::transform(dev_vec.begin(), dev_vec.end(), dev_vec.begin(), factor * thrust::placeholders::_1);
  }
}

template<typename xray::KernelPrecision PRECISION>
xray::DPartialGPU<PRECISION>::DPartialGPU(int n_hkl, const int* hkl, const double* mss4, std::complex<double>* f_calc,
                               const double* abs_f_calc, int n_atom, const double* atom_b_factor, const double* atom_occupancy,
                               const int* atom_scatter_type, int n_scatter_types, const double* atomic_scatter_factor)
  : xray::DPartial(n_hkl, hkl, mss4, f_calc, abs_f_calc, n_atom, atom_b_factor, atom_occupancy, atom_scatter_type,
                   n_scatter_types, atomic_scatter_factor) {

  m_dev_frac_by_2_pi = thrust::device_vector<FloatType>(n_atom * 3);
  m_dev_d_target_d_abs_f_calc_by_2_pi = thrust::device_vector<FloatType>(n_hkl);
  m_dev_atomic_scatter_factor = thrust::device_vector<FloatType>(m_atomic_scatter_factor,
                                                              m_atomic_scatter_factor + m_n_scatter_types * n_hkl);
  m_dev_abs_f_calc = thrust::device_vector<FloatType>(n_hkl);
  m_dev_f_calc_phase = thrust::device_vector<FloatType>(n_hkl);
  m_dev_f_scale = thrust::device_vector<FloatType>(n_hkl);
  m_dev_hkl = thrust::device_vector<int>(m_hkl, m_hkl + n_hkl * 3);
  m_dev_mss4 = thrust::device_vector<FloatType>(m_mss4, m_mss4 + n_hkl);
  m_dev_atom_b_factor = thrust::device_vector<FloatType>(m_atom_b_factor, m_atom_b_factor + n_atom);
  m_dev_atom_occupancy = thrust::device_vector<FloatType>(m_atom_occupancy, m_atom_occupancy + n_atom);
  m_dev_atom_scatter_type = thrust::device_vector<int>(m_atom_scatter_type, m_atom_scatter_type + n_atom);
  m_dev_d_target_d_frac = thrust::device_vector<FloatType>(n_atom * 3);
  m_f_calc_phase = std::vector<FloatType>(n_hkl);
}

template<typename xray::KernelPrecision PRECISION>
void xray::DPartialGPU<PRECISION>::calc_d_target_d_frac(
  int n_atom,
  const double* frac,
  int n_hkl,
  const double* f_scale,
  const double* d_target_d_abs_f_calc,
  double* d_target_d_frac) {

  assert(n_atom == m_n_atom);
  assert(n_hkl == m_n_hkl);

  thrust::copy(frac, frac + n_atom * 3, m_dev_frac_by_2_pi.begin());
  thrust::copy(m_abs_f_calc, m_abs_f_calc + n_hkl, m_dev_abs_f_calc.begin());
  multiply_by(static_cast<FloatType>(2 * M_PI), m_dev_frac_by_2_pi);
  thrust::copy(d_target_d_abs_f_calc, d_target_d_abs_f_calc + n_hkl, m_dev_d_target_d_abs_f_calc_by_2_pi.begin());
  thrust::copy(f_scale, f_scale + n_hkl, m_dev_f_scale.begin());
  multiply_by(static_cast<FloatType>(2 * M_PI), m_dev_d_target_d_abs_f_calc_by_2_pi);

  for (int i = 0; i < n_hkl; i++) {
    m_f_calc_phase[i] = arg(m_f_calc[i]);
  }

  thrust::copy(m_f_calc_phase.begin(), m_f_calc_phase.end(), m_dev_f_calc_phase.begin());

  dim3 numBlocks(n_atom);
  const int block_size = 64;
  dim3 threadsPerBlock(block_size);

  calc_d_target_d_frac_kernel<block_size>
  <<<numBlocks, threadsPerBlock>>>(
    n_atom,
    m_dev_frac_by_2_pi.data().get(),
    m_dev_atom_b_factor.data().get(),
    m_dev_atom_occupancy.data().get(),
    m_n_hkl,
    m_dev_hkl.data().get(),
    m_dev_mss4.data().get(),
    m_dev_atomic_scatter_factor.data().get(),
    m_dev_atom_scatter_type.data().get(),
    m_dev_f_calc_phase.data().get(),
    m_dev_f_scale.data().get(),
    m_dev_abs_f_calc.data().get(),
    m_dev_d_target_d_abs_f_calc_by_2_pi.data().get(),
    m_dev_d_target_d_frac.data().get()
  );

  thrust::copy(m_dev_d_target_d_frac.begin(), m_dev_d_target_d_frac.end(), d_target_d_frac);
}

template class xray::DPartialGPU<xray::KernelPrecision::Single>;
template class xray::DPartialGPU<xray::KernelPrecision::Double>;
