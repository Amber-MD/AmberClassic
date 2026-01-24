module xray_target_vector_least_squares_impl_cpu_module

    use xray_contracts_module
    use xray_pure_utils, only : real_kind
    use xray_interface2_data_module, only : n_hkl, Fobs, n_work
    use xray_target_vector_least_squares_data_module

    implicit none

    private

    ! Public module interface
    public :: init
    public :: calc_partial_d_target_d_Fcalc
    public :: finalize

contains

    subroutine init(abs_Fobs)
        implicit none
        real(real_kind), intent(in) :: abs_Fobs(:)
        integer i

        allocate(derivc(n_hkl))
        norm_scale = 1 / sum(abs(Fobs(:n_work)) ** 2)
        write(6,'(a,e14.7)') '| vls: setting norm_scale = ', norm_scale

    end subroutine init

    subroutine finalize()
    end subroutine finalize

    ! This routine computes the force gradient on Fcalc as a harmonic
    ! restraint on the vector (complex) difference between Fcalc and
    ! Fobs

    ! derivc is treated as a complex variable, but this is not really
    ! complex analysis!  Rather, derivc%re is the derivative of the target
    ! with respect to Fcalc%re, and derivc%im is the derivative of the 
    ! target with respect to Fcalc%im.

    subroutine calc_partial_d_target_d_Fcalc(Fcalc, derivc, xray_energy)
        implicit none
        complex(real_kind), intent(in) :: Fcalc(:)
        complex(real_kind), intent(out) :: derivc(:)
        real(real_kind), intent(out) :: xray_energy
        complex(real_kind) :: vecdif(size(Fcalc))

        integer :: i

        vecdif(:) = Fobs(:) - Fcalc(:)
        xray_energy = norm_scale * sum(vecdif(:n_work) * conjg(vecdif(:n_work)))

        derivc(:n_work) = - norm_scale * 2 * vecdif(:n_work)
        derivc(n_work + 1:) = 0 ! no force for things unselected here


    end subroutine calc_partial_d_target_d_Fcalc

end module xray_target_vector_least_squares_impl_cpu_module
