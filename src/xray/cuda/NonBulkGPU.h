#ifndef AMBER_NON_BULK_GPU_H
#define AMBER_NON_BULK_GPU_H

#include "NonBulk.h"
#include <thrust/device_vector.h>

namespace xray {

  enum class NonBulkKernelVersion {
    ManualCaching,
    StraightForward,
  };

  template<NonBulkKernelVersion KERNEL_VERSION, KernelPrecision PRECISION>
  class NonBulkGPU : public NonBulk {
    using FloatType = typename KernelConfig<PRECISION>::FloatType;
  public:
    NonBulkGPU(int n_hkl,
               const int* hkl,
               complex_double* f_non_bulk,
               const double* mSS4,
               int n_atom,
               int nb,
               const double* b_factor,
               const double* occupancy,
               int n_scatter_types,
               const int* scatter_type_index,
               const double* atomic_scatter_factor);

    ~NonBulkGPU() override = default;

    void calc_f_non_bulk(
      int n_atoms,
      const double* frac_xyz
    ) override;

  private:
    thrust::device_vector<FloatType> m_dev_frac_xyz;
    thrust::device_vector<FloatType> m_dev_b_factor;
    thrust::device_vector<FloatType> m_dev_occupancy;
    thrust::device_vector<int> m_dev_hkl;
    thrust::device_vector<FloatType> m_dev_mSS4;
    thrust::device_vector<FloatType> m_dev_atomic_scatter_factor;

    thrust::device_vector<int> m_dev_scatter_type_index;
    thrust::device_vector<thrust::complex<FloatType>> m_dev_f_non_bulk;

  };
}
#endif //AMBER_NON_BULK_GPU_H
