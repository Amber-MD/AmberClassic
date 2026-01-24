#ifndef GpuPortingDataStructures
#define GpuPortingDataStructures

// Definitions relevant to CUDA kernels and fixed precision arithmetic
#define FPSCALErad      (float)16777216.0
#define FPSCALEfrc      (float)524288.0
#define FPSCALEnrgHP    (float)524288.0
#define FPSCALEnrgMP    (float)131072.0
#define FPSCALEnrgLP    (float)32768.0
#define FPSCALEcn       (float)1073741824.0
#define FPSCALEicn      (float)9.3132257461547852e-10
#define FPSCALEcnv      (float)268435456.0
#define FPSCALEicnv     (float)3.7252902984619140625e-09
#define RCSCALE         (float)512.0
#define HP2LP_BITS      4
#define HP2MP_BITS      2
#define MP2LP_BITS      2
#if defined(AMBER_PLATFORM_AMD)
#define QUARTGRID        16
#define HALFGRID         32
#define GRID             64
#define GRIDx2           128
#define GRIDx3           192
#define GRIDx4           256
#define GRIDx5           320
#define GRIDx6           384
#define HALFGRID_MASK    31
#define GRID_BITS_MASK   63
#define GRID_BITS_DIAG   5
#define GRID_BITS        6
#else
#define QUARTGRID        8
#define HALFGRID         16
#define GRID             32
#define GRIDx2           64
#define GRIDx3           96
#define GRIDx4           128
#define GRIDx5           160
#define GRIDx6           192
#define HALFGRID_MASK    15
#define GRID_BITS_MASK   31
#define GRID_BITS_DIAG   4
#define GRID_BITS        5
#endif

// Taylor series coefficients for GB
#define TA  (float)0.33333333333333333333
#define TB  (float)0.4
#define TC  (float)0.42857142857142857143
#define TD  (float)0.44444444444444444444
#define TDD (float)0.45454545454545454545
#define TE  (float)1.33333333333333333333
#define TF  (float)2.4
#define TG  (float)3.42857142857142857143
#define TH  (float)4.44444444444444444444
#define THH (float)5.45454545454545454545

#  ifndef CUDA
// If CUDA is not defined, then some of the CUDA POD vector types will need to
// be delineated in order to run the program in CPU code.
struct IntTwo {
  int x;
  int y;
};
typedef struct IntTwo int2;

struct IntFour {
  int x;
  int y;
  int z;
  int w;
};
typedef struct IntFour int4;

typedef unsigned int uint;

struct UIntTwo {
  unsigned int x;
  unsigned int y;
};
typedef struct UIntTwo uint2;

struct UIntFour {
  unsigned int x;
  unsigned int y;
  unsigned int z;
  unsigned int w;
};
typedef struct UIntFour uint4;

struct FloatTwo {
  float x;
  float y;
};
typedef struct FloatTwo float2;

struct FloatFour {
  float x;
  float y;
  float z;
  float w;
};
typedef struct FloatFour float4;

struct DoubleTwo {
  double x;
  double y;
};
typedef struct DoubleTwo double2;

struct DoubleFour {
  double x;
  double y;
  double z;
  double w;
};
typedef struct DoubleFour double4;
#  else
#    if defined(AMBER_PLATFORM_AMD)
       // Must be replaced as soon as ROCm is updated
#      include <hip/hip_runtime_api.h>
#      include "HipDefinitions.h"
#    else
#      include <cuda.h>
#      include <cuda_runtime_api.h>
#    endif
#  endif

struct GpuMirroredInt {
  int len;
  int IsPinned;
  int* HostData;
  int* DevcData;
};
typedef struct GpuMirroredInt gpuInt;

struct GpuMirroredInt2 {
  int len;
  int IsPinned;
  int2* HostData;
  int2* DevcData;
};
typedef struct GpuMirroredInt2 gpuInt2;

struct GpuMirroredInt4 {
  int len;
  int IsPinned;
  int4* HostData;
  int4* DevcData;
};
typedef struct GpuMirroredInt4 gpuInt4;

struct GpuMirroredUInt {
  int len;
  int IsPinned;
  uint* HostData;
  uint* DevcData;
};
typedef struct GpuMirroredUInt gpuUInt;

struct GpuMirroredUInt2 {
  int len;
  int IsPinned;
  uint2* HostData;
  uint2* DevcData;
};
typedef struct GpuMirroredUInt2 gpuUInt2;

struct GpuMirroredUInt4 {
  int len;
  int IsPinned;
  uint4* HostData;
  uint4* DevcData;
};
typedef struct GpuMirroredUInt4 gpuUInt4;

struct GpuMirroredFloat {
  int len;
  int IsPinned;
  float* HostData;
  float* DevcData;
};
typedef struct GpuMirroredFloat gpuFloat;

struct GpuMirroredFloat2 {
  int len;
  int IsPinned;
  float2* HostData;
  float2* DevcData;
};
typedef struct GpuMirroredFloat2 gpuFloat2;

struct GpuMirroredFloat4 {
  int len;
  int IsPinned;
  float4* HostData;
  float4* DevcData;
};
typedef struct GpuMirroredFloat4 gpuFloat4;

struct GpuMirroredDouble {
  int len;
  int IsPinned;
  double* HostData;
  double* DevcData;
};
typedef struct GpuMirroredDouble gpuDouble;

struct GpuMirroredDouble2 {
  int len;
  int IsPinned;
  double2* HostData;
  double2* DevcData;
};
typedef struct GpuMirroredDouble2 gpuDouble2;

struct GpuMirroredDouble4 {
  int len;
  int IsPinned;
  double4* HostData;
  double4* DevcData;
};
typedef struct GpuMirroredDouble4 gpuDouble4;

#endif
