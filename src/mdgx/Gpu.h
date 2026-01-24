#ifndef GpuMirroringFunctions
#define GpuMirroringFunctions

#include "GpuDS.h"

gpuInt CreateGpuInt(int len, int pin);

void DestroyGpuInt(gpuInt *G);

gpuInt CopyGpuInt(gpuInt *G, int copySrc);

void UploadGpuInt(gpuInt *G, int* other);

void DownloadGpuInt(gpuInt *G, int* other);

gpuInt2 CreateGpuInt2(int len, int pin);

void DestroyGpuInt2(gpuInt2 *G);

gpuInt2 CopyGpuInt2(gpuInt2 *G, int copySrc);

void UploadGpuInt2(gpuInt2 *G, int2* other);

void DownloadGpuInt2(gpuInt2 *G, int2* other);

gpuInt4 CreateGpuInt4(int len, int pin);

void DestroyGpuInt4(gpuInt4 *G);

gpuInt4 CopyGpuInt4(gpuInt4 *G, int copySrc);

void UploadGpuInt4(gpuInt4 *G, int4* other);

void DownloadGpuInt4(gpuInt4 *G, int4* other);

gpuUInt CreateGpuUInt(int len, int pin);

void DestroyGpuUInt(gpuUInt *G);

gpuUInt CopyGpuUInt(gpuUInt *G, int copySrc);

void UploadGpuUInt(gpuUInt *G, uint* other);

void DownloadGpuUInt(gpuUInt *G, uint* other);

gpuUInt2 CreateGpuUInt2(int len, int pin);

void DestroyGpuUInt2(gpuUInt2 *G);

gpuUInt2 CopyGpuUInt2(gpuUInt2 *G, int copySrc);

void UploadGpuUInt2(gpuUInt2 *G, uint2* other);

void DownloadGpuUInt2(gpuUInt2 *G, uint2* other);

gpuUInt4 CreateGpuUInt4(int len, int pin);

void DestroyGpuUInt4(gpuUInt4 *G);

gpuUInt4 CopyGpuUInt4(gpuUInt4 *G, int copySrc);

void UploadGpuUInt4(gpuUInt4 *G, uint4* other);

void DownloadGpuUInt4(gpuUInt4 *G, uint4* other);

gpuFloat CreateGpuFloat(int len, int pin);

void DestroyGpuFloat(gpuFloat *G);

gpuFloat CopyGpuFloat(gpuFloat *G, int copySrc);

void UploadGpuFloat(gpuFloat *G, float* other);

void DownloadGpuFloat(gpuFloat *G, float* other);

gpuFloat2 CreateGpuFloat2(int len, int pin);

void DestroyGpuFloat2(gpuFloat2 *G);

gpuFloat2 CopyGpuFloat2(gpuFloat2 *G, int copySrc);

void UploadGpuFloat2(gpuFloat2 *G, float2* other);

void DownloadGpuFloat2(gpuFloat2 *G, float2* other);

gpuFloat4 CreateGpuFloat4(int len, int pin);

void DestroyGpuFloat4(gpuFloat4 *G);

gpuFloat4 CopyGpuFloat4(gpuFloat4 *G, int copySrc);

void UploadGpuFloat4(gpuFloat4 *G, float4* other);

void DownloadGpuFloat4(gpuFloat4 *G, float4* other);

gpuDouble CreateGpuDouble(int len, int pin);

void DestroyGpuDouble(gpuDouble *G);

gpuDouble CopyGpuDouble(gpuDouble *G, int copySrc);

void UploadGpuDouble(gpuDouble *G, double* other);

void DownloadGpuDouble(gpuDouble *G, double* other);

gpuDouble2 CreateGpuDouble2(int len, int pin);

void DestroyGpuDouble2(gpuDouble2 *G);

gpuDouble2 CopyGpuDouble2(gpuDouble2 *G, int copySrc);

void UploadGpuDouble2(gpuDouble2 *G, double2* other);

void DownloadGpuDouble2(gpuDouble2 *G, double2* other);

gpuDouble4 CreateGpuDouble4(int len, int pin);

void DestroyGpuDouble4(gpuDouble4 *G);

gpuDouble4 CopyGpuDouble4(gpuDouble4 *G, int copySrc);

void UploadGpuDouble4(gpuDouble4 *G, double4* other);

void DownloadGpuDouble4(gpuDouble4 *G, double4* other);

#endif
