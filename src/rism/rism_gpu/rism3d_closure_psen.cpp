#include <iostream>
#include "rism3d_closure_psen.hpp"
#include <sstream>
#include <iostream>
using namespace std;

namespace rism3d_c {

    rism3d_closure_psen :: rism3d_closure_psen(){
        cout << "Default object from rism3d_closure_psen is being created!" << endl;
    }

    rism3d_closure_psen :: rism3d_closure_psen(rism3d_solvent_class *solventclass, rism3d_solute_class *soluteclass, rism3d_grid *grid, rism3d_potential *pot, rism3d_safemem *memalloc){
        cout << "Object from rism3d_closure_psen is being created!" << endl;
        solv = solventclass;
    }

    rism3d_closure_psen :: ~rism3d_closure_psen(){
        cout << "Object from rism3d_closure_psen is being deleted!" << endl;
    }

    void rism3d_closure_psen :: guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size){
        cout << "from PSE-n class: Hello" << endl;
    }  

}