//-----------------------------------------------------------------------------
// Yes, mdgx is not C++.  This include will replicate the same functions for
// many data types, in what should be a template class in C++.  The following
// macros must be #define'd when including this code:
//
// MIRTYPE:     the mirror memoy struct type (e.g. gpuInt)
// CONSTRUCTOR: the name for the constructor function
// ELMTYPE:     the type of each element in the array to allocate (e.g. int)
// DESTRUCTOR:  the name for the destructor function
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// CONSTRUCTOR: create an array mirrored on both host and device memory.  This
//              is typically going to be called CreateMirrored<Data Type>().
//
// Arguments:
//   len:      the length of array to allocate
//   pin:      flag to have the memory pinned (non-pageable on the host side
//             for optimal transfer speed ot the device)
//-----------------------------------------------------------------------------
MIRTYPE CONSTRUCTOR(int len, int pin)
{
  MIRTYPE G;

  G.len = len;
  G.IsPinned = pin;
#ifdef CUDA
  len = ((len + 127) / 128) * 128;
  if (pin == 1) {
    cudaHostAlloc((void **)&G.HostData, len * sizeof(PODTYPE),
                  cudaHostAllocMapped);
  }
  else {
    G.HostData = (PODTYPE*)malloc(len * sizeof(PODTYPE));
  }
  cudaMalloc((void **)&G.DevcData, len * sizeof(PODTYPE));
  memset(G.HostData, 0, len * sizeof(PODTYPE));
  cudaMemset((void *)G.DevcData, 0, len * sizeof(PODTYPE));
#else
  G.HostData = (PODTYPE*)calloc(len, sizeof(PODTYPE));
#endif

  return G;
}

//-----------------------------------------------------------------------------
// DESTRUCTOR: free an array mirrored on both host and device memory.  This
//             is typically going to be called DestroyMirrored<Data Type>().
//
// Arguments:
//   G:    the gpu<data type> struct to free
//-----------------------------------------------------------------------------
void DESTRUCTOR(MIRTYPE *G)
{
#ifdef CUDA
  if (G->IsPinned == 1) {
    cudaFreeHost(G->HostData);
  }
  else {
    free(G->HostData);
  }
  cudaFree(G->DevcData);
#else
  free(G->HostData);
#endif
}

//-----------------------------------------------------------------------------
// COPYRTN: function to copy a struct containing GPU mirrored memory.
//
// Arguments:
//   G:        the original gpu<data type> struct to free
//   copySrc:  set to 0 to copy from host memory, 1 to copy from device memory
//-----------------------------------------------------------------------------
MIRTYPE COPYRTN(MIRTYPE *G, int copySrc)
{
  int i;
  MIRTYPE newG;
#ifdef CUDA
  cudaError_t status; 
#endif
  
  newG = CONSTRUCTOR(G->len, G->IsPinned);
  if (copySrc == 0) {
    for (i = 0; i < G->len; i++) {
      newG.HostData[i] = G->HostData[i];
    }
  }
  else if (copySrc == 1) {
#ifdef CUDA
    status = cudaMemcpy(G->DevcData, newG.HostData, G->len * sizeof(MIRTYPE),
                        cudaMemcpyHostToDevice);
    if (status != cudaSuccess) {
      printf("mdgx >> Error.  Copy for mirrored array of length %d failed.\n",
	     G->len);
    }
#else
    return newG;
#endif
  }
  
  return newG;
}

//-----------------------------------------------------------------------------
// The following functions are only supported if CUDA is in the compilation.
//-----------------------------------------------------------------------------
#ifdef CUDA

//-----------------------------------------------------------------------------
// UPLOADER: function to upload memory in a gpu<data type> struct to the GPU
//           from host RAM.
//
// Arguments:
//   G:  the struct to upload
//-----------------------------------------------------------------------------
void UPLOADER(MIRTYPE *G, PODTYPE *other)
{
#ifdef CUDA
  cudaError_t status;
#endif
  
  if (other != NULL) {
    status = cudaMemcpy(G->DevcData, other, G->len * sizeof(PODTYPE),
                        cudaMemcpyHostToDevice);
  }
  else {
    status = cudaMemcpy(G->DevcData, G->HostData, G->len * sizeof(PODTYPE),
                        cudaMemcpyHostToDevice);
  }
}

//-----------------------------------------------------------------------------
// DOWNLOADER: function to download memory in a gpu<data type> struct from the
//             GPU back to host RAM.
//
// Arguments:
//   G:  the struct to upload
//-----------------------------------------------------------------------------
void DOWNLOADER(MIRTYPE *G, PODTYPE *other)
{
  cudaError_t status;

  if (other != NULL) {
    status = cudaMemcpy(other, G->DevcData, G->len * sizeof(PODTYPE),
                        cudaMemcpyDeviceToHost);
  }
  else {
    status = cudaMemcpy(G->HostData, G->DevcData, G->len * sizeof(PODTYPE),
                        cudaMemcpyDeviceToHost);
  }
}
#endif
