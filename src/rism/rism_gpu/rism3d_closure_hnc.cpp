#include <iostream>
#include "rism3d_closure_hnc.hpp"
#include <sstream>
#include <iostream>
using namespace std;

namespace rism3d_c {

    rism3d_closure_hnc :: rism3d_closure_hnc(){
        cout << "Default object from rism3d_closure_hnc is being created!" << endl;
    }

    rism3d_closure_hnc :: rism3d_closure_hnc(rism3d_solvent_class *solventclass, rism3d_solute_class *soluteclass, rism3d_grid *grid, rism3d_potential *pot, rism3d_safemem *memalloc){
        cout << "Object from rism3d_closure_hnc is being created!" << endl;
        solv = solventclass;
    }

    rism3d_closure_hnc :: ~rism3d_closure_hnc(){
        cout << "Object from rism3d_closure_hnc is being deleted!" << endl;
    }

    void rism3d_closure_hnc :: guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size){
        cout << "from HNC class: Hello" << endl;
    }  


}