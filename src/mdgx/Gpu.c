#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Gpu.h"
#include "mdgxVector.h"

//-----------------------------------------------------------------------------
// Functions for allocating and deallocating mirrored memory on the host and
// device.  These draw from the GpuBuffer template class in pmemd, but
// simplify the mechanics a tad.
//-----------------------------------------------------------------------------
#define MIRTYPE gpuInt
#define PODTYPE int
#define CONSTRUCTOR CreateGpuInt
#define DESTRUCTOR DestroyGpuInt
#define COPYRTN CopyGpuInt
#define UPLOADER UploadGpuInt
#define DOWNLOADER DownloadGpuInt
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuInt2
#define PODTYPE int2
#define CONSTRUCTOR CreateGpuInt2
#define DESTRUCTOR DestroyGpuInt2
#define COPYRTN CopyGpuInt2
#define UPLOADER UploadGpuInt2
#define DOWNLOADER DownloadGpuInt2
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuInt4
#define PODTYPE int4
#define CONSTRUCTOR CreateGpuInt4
#define DESTRUCTOR DestroyGpuInt4
#define COPYRTN CopyGpuInt4
#define UPLOADER UploadGpuInt4
#define DOWNLOADER DownloadGpuInt4
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuUInt
#define PODTYPE uint
#define CONSTRUCTOR CreateGpuUInt
#define DESTRUCTOR DestroyGpuUInt
#define COPYRTN CopyGpuUInt
#define UPLOADER UploadGpuUInt
#define DOWNLOADER DownloadGpuUInt
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuUInt2
#define PODTYPE uint2
#define CONSTRUCTOR CreateGpuUInt2
#define DESTRUCTOR DestroyGpuUInt2
#define COPYRTN CopyGpuUInt2
#define UPLOADER UploadGpuUInt2
#define DOWNLOADER DownloadGpuUInt2
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuUInt4
#define PODTYPE uint4
#define CONSTRUCTOR CreateGpuUInt4
#define DESTRUCTOR DestroyGpuUInt4
#define COPYRTN CopyGpuUInt4
#define UPLOADER UploadGpuUInt4
#define DOWNLOADER DownloadGpuUInt4
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuFloat
#define PODTYPE float
#define CONSTRUCTOR CreateGpuFloat
#define DESTRUCTOR DestroyGpuFloat
#define COPYRTN CopyGpuFloat
#define UPLOADER UploadGpuFloat
#define DOWNLOADER DownloadGpuFloat
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuFloat2
#define PODTYPE float2
#define CONSTRUCTOR CreateGpuFloat2
#define DESTRUCTOR DestroyGpuFloat2
#define COPYRTN CopyGpuFloat2
#define UPLOADER UploadGpuFloat2
#define DOWNLOADER DownloadGpuFloat2
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuFloat4
#define PODTYPE float4
#define CONSTRUCTOR CreateGpuFloat4
#define DESTRUCTOR DestroyGpuFloat4
#define COPYRTN CopyGpuFloat4
#define UPLOADER UploadGpuFloat4
#define DOWNLOADER DownloadGpuFloat4
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuDouble
#define PODTYPE double
#define CONSTRUCTOR CreateGpuDouble
#define DESTRUCTOR DestroyGpuDouble
#define COPYRTN CopyGpuDouble
#define UPLOADER UploadGpuDouble
#define DOWNLOADER DownloadGpuDouble
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuDouble2
#define PODTYPE double2
#define CONSTRUCTOR CreateGpuDouble2
#define DESTRUCTOR DestroyGpuDouble2
#define COPYRTN CopyGpuDouble2
#define UPLOADER UploadGpuDouble2
#define DOWNLOADER DownloadGpuDouble2
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER

#define MIRTYPE gpuDouble4
#define PODTYPE double4
#define CONSTRUCTOR CreateGpuDouble4
#define DESTRUCTOR DestroyGpuDouble4
#define COPYRTN CopyGpuDouble4
#define UPLOADER UploadGpuDouble4
#define DOWNLOADER DownloadGpuDouble4
#include "GpuMirrors.h"
#undef MIRTYPE
#undef PODTYPE
#undef CONSTRUCTOR
#undef DESTRUCTOR
#undef COPYRTN
#undef UPLOADER
#undef DOWNLOADER
