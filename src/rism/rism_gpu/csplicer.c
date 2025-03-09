// splicer begin namespace.rism3d_c.class.rism3d.CXX_declarations
#include "wrapRISM3D_rism3d_c.h"
// splicer end namespace.rism3d_c.class.rism3d.CXX_declarations

///////////////////////////////////////////////////////////////////////

// splicer begin namespace.rism3d_c.class.rism3d.C_declarations
int RIS_rism3d_c_rism3d_get_numatoms(RIS_rism3d_c_rism3d * self);

int RIS_rism3d_c_rism3d_get_numatomtypes(RIS_rism3d_c_rism3d * self);
// splicer end namespace.rism3d_c.class.rism3d.C_declarations

// splicer begin namespace.rism3d_c.class.rism3d.C_definitions
int RIS_rism3d_c_rism3d_get_numatoms(RIS_rism3d_c_rism3d * self)
{
    rism3d_c::rism3d *SH_this =
        static_cast<rism3d_c::rism3d *>(self->addr);
    int SHC_rv = SH_this->get_numatoms();
    return SHC_rv;
}

int RIS_rism3d_c_rism3d_get_numatomtypes(RIS_rism3d_c_rism3d * self)
{
    rism3d_c::rism3d *SH_this = static_cast<rism3d_c::rism3d *>
        (self->addr);
    int SHC_rv = SH_this->get_numatomtypes();
    return SHC_rv;
}
// splicer end namespace.rism3d_c.class.rism3d.C_definitions