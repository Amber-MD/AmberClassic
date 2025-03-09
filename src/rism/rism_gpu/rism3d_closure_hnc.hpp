#ifndef RISM3DCLOSUREHNC_HPP
#define RISM3DCLOSUREHNC_HPP

#include <iostream>
#include "rism3d_closure.hpp"
using namespace std;

namespace rism3d_c {

    class rism3d_closure_hnc : public rism3d_closure {            
        public:
            rism3d_closure_hnc();
            rism3d_closure_hnc(rism3d_solvent_class *solventclass, rism3d_solute_class *soluteclass, 
                           rism3d_grid *grid, rism3d_potential *pot, rism3d_safemem *memalloc);
            ~rism3d_closure_hnc();
            
            void guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size) override;

    };

}

#endif // RISM3DCLOSUREHNC_HPP