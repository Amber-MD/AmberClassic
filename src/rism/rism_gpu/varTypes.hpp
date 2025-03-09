#ifndef VARTYPES_HPP

#if defined(RISMCUDA_DOUBLE)
typedef double GPUtype;
#else
typedef float GPUtype;
#endif

#endif // VARTYPES_HPP