! splicer begin namespace.rism3d_c.module_use
use rism3d_solute_c
use rism3d_solvent_c
! splicer end namespace.rism3d_c.module_use

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! splicer begin namespace.rism3d_c.additional_declarations
interface
    
    pure function c_rism3d_get_numatoms(self) &
            result(SHT_rv) &
            bind(C, name="RIS_rism3d_c_rism3d_get_numatoms")
        use iso_c_binding, only : C_INT
        import :: RIS_SHROUD_capsule_data
        implicit none
        type(RIS_SHROUD_capsule_data), intent(IN) :: self
        integer(C_INT) :: SHT_rv
    end function c_rism3d_get_numatoms

    pure function c_rism3d_get_numatomtypes(self) &
        result(SHT_rv) &
        bind(C, name="RIS_rism3d_c_rism3d_get_numatomtypes")
        use iso_c_binding, only : C_INT
        import :: RIS_SHROUD_capsule_data
        implicit none
        type(RIS_SHROUD_capsule_data), intent(IN) :: self
        integer(C_INT) :: SHT_rv
    end function c_rism3d_get_numatomtypes

end interface

interface  rism3d_map_site_to_site
    module procedure rism3d_map_site_to_site_flat, rism3d_map_site_to_site_3_D
end interface rism3d_map_site_to_site
! splicer end namespace.rism3d_c.additional_declarations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! splicer begin namespace.rism3d_c.class.rism3d.additional_functions
pure function rism3d_get_numatoms(obj) &
        result(SHT_rv)
    use iso_c_binding, only : C_INT
    class(rism3d), intent(in) :: obj
    integer(C_INT) :: SHT_rv
    SHT_rv = c_rism3d_get_numatoms(obj%cxxmem)
end function rism3d_get_numatoms

pure function rism3d_get_numatomtypes(obj) &
        result(SHT_rv)
    use iso_c_binding, only : C_INT
    class(rism3d), intent(in) :: obj
    integer(C_INT) :: SHT_rv
    SHT_rv = c_rism3d_get_numatomtypes(obj%cxxmem)
end function rism3d_get_numatomtypes

! splicer end namespace.rism3d_c.class.rism3d.additional_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! splicer begin namespace.rism3d_c.class.rism3d.type_bound_procedure_part

! splicer end namespace.rism3d_c.class.rism3d.type_bound_procedure_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! splicer begin namespace.rism3d_c.additional_functions
subroutine cpy_solute(this, solute)
    type(solute_cpp) :: this
    type(rism3d_solute) :: solute
    REAL, target :: COM(3)
    REAL, target :: trans(3)

    COM = solute%centerOfMass
    trans = solute%translation
    
    this%numAtoms = solute%numAtoms
    this%charged = solute%charged
    this%mass = C_LOC(solute%mass)
    this%charge = C_LOC(solute%charge)
    this%ljSigma = C_LOC(solute%ljSigma)
    this%ljEpsilon = C_LOC(solute%ljEpsilon)
    this%centerOfMass = C_LOC(COM)
    this%translation = C_LOC(trans)
    this%origCharge = C_LOC(solute%origCharge)
    this%position = C_LOC(solute%position)
    this%totalCharge = solute%totalCharge
end subroutine cpy_solute

subroutine cpy_solvent(this, solvent)
    type(solvent_cpp) :: this
    type(rism3d_solvent) :: solvent
    
    this%temperature = solvent%temperature
    this%dielconst = solvent%dielconst
    this%xappa = solvent%xappa
    this%xikt = solvent%xikt
    this%smear = solvent%smear
    this%xikt_dT = solvent%xikt_dT
    this%numAtomTypes = solvent%numAtomTypes
    this%numMolecules = solvent%numMolecules
    this%numRDFpoints = solvent%numRDFpoints
    this%atomMultiplicity = C_LOC(solvent%atomMultiplicity)
    this%numAtoms = C_LOC(solvent%numAtoms)
    this%gridSpacingR = solvent%gridSpacingR
    this%gridSpacingK = solvent%gridSpacingK
    this%waveNumbers = C_LOC(solvent%waveNumbers)
    this%xvv = C_LOC(solvent%xvv)
    this%xvv_dT = C_LOC(solvent%xvv_dT)
    this%charge = C_LOC(solvent%charge)
    this%charge_sp = C_LOC(solvent%charge_sp)
    this%density = C_LOC(solvent%density)
    this%density_sp = C_LOC(solvent%density_sp)
    this%ljSigma = C_LOC(solvent%ljSigma)
    this%ljEpsilon = C_LOC(solvent%ljEpsilon)
    this%coord = C_LOC(solvent%coord)
    this%background_correction = C_LOC(solvent%background_correction)
    this%delhv0 = C_LOC(solvent%delhv0)
    this%delhv0_dT = C_LOC(solvent%delhv0_dT)
    this%ionic = solvent%ionic
    this%xvv_version = solvent%xvv_version
end subroutine cpy_solvent
! splicer end namespace.rism3d_c.additional_functions