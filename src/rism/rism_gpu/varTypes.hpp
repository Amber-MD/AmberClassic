#ifndef VARTYPES_HPP
#define VARTYPES_HPP

#if defined(RISMCUDA_DOUBLE)

    #define GPUtype double
    #define GPUReduceAccumType double
    #define GPUMDIISType double
    #define GPUPotAccumType double

#else

    #define GPUtype float
    // #define GPUReduceAccumType double
    #define GPUReduceAccumType float
    #define GPUMDIISType double
    // #define GPUMDIISType float
    #define GPUPotAccumType double
    // #define GPUPotAccumType float
#endif

#endif // VARTYPES_HPP