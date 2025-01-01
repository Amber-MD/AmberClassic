//-----------------------------------------------------------------------------
// Useful macros for kernel operations on the GPU
//-----------------------------------------------------------------------------
#ifndef KernelMacroFramework
#define KernelMacroFramework

#if defined(AMBER_PLATFORM_AMD) || (defined(CUDA_VERSION) && (CUDA_VERSION < 9000))
#  define __shfl_down_sync(a, b, c, d)  __shfl_down(b, c, d)
#  define __shfl_xor_sync(a, b, c, d)   __shfl_xor(b, c, d)
#  define __shfl_up_sync(a, b, c, d)    __shfl_up(b, c, d)
#  define __shfl_sync(a, b, c, d)       __shfl(b, c, d)
#  define __ballot_sync(a, b)           __ballot(b)
#  define __syncwarp(a)
#endif

#  include "GpuDS.h"

#if defined(AMBER_PLATFORM_AMD)
__device__ inline float fexp(float x)
{
    return __expf(x);
}

__device__ inline float flog(float x)
{
    return __logf(x);
}

__device__ inline float fsqrt(const float x)
{
    // using hardware instruction causes some tests to fail by a small margin
    // return __fsqrt_rn(x);
    // so we'll use the babylonian method
    union
    {
        int i;
        float x;
    } u;
    u.x = x;
    u.i = (1 << 29) + (u.i >> 1) - (1 << 22);
    u.x = u.x + x / u.x;
    u.x = 0.25f * u.x + x / u.x;
    return u.x;
}

__device__ inline float frsqrt(const float x)
{
    return __frsqrt_rn(x);
}

#define exp fexp
#define sqrt fsqrt
#define rsqrt frsqrt
#define log flog
#endif

//-----------------------------------------------------------------------------
// WarpREDUCE: warp reduction to sum all values into the first lane
//-----------------------------------------------------------------------------
#define WarpREDUCE(var) \
{ \
  for(unsigned int offset = GRID >> 1; offset > 0; offset >>= 1) \
  { \
      var += __shfl_down_sync(0xffffffff, var, offset, GRID); \
  } \
}

//-----------------------------------------------------------------------------
// DvcCROSSPf: compute the cross product vA x vB = vC
//-----------------------------------------------------------------------------
#define DvcCROSSPf(vA, vB, vC) \
{ \
  vC[0] = vA[1]*vB[2] - vA[2]*vB[1]; \
  vC[1] = vA[2]*vB[0] - vA[0]*vB[2]; \
  vC[2] = vA[0]*vB[1] - vA[1]*vB[0]; \
}

#define ISNAN(x) ((__float_as_uint(x) & 0x7f000000) == 0x7f000000 &&	\
		  (__float_as_uint(x) & 0xffffff) > 0)

#endif
