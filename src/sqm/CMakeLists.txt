project(sqm)

set(EXECUTABLE_SOURCES sqm.F90 
	xmin.F90
	amopen.F90 mexit.F90 assert.F90 
	timer_dummy.F90 nmlsrc.F90 
	qm_print_info.F90 
	qm2_energy.F90 
	qm2_read_nm_and_alloc.F90 
	qm2_scf.F90
	qm2_dftb_scf.F90)

	
set(LIBRARY_SOURCES qm_print_info.F90 
	qm2_energy.F90 
	qm2_read_nm_and_alloc.F90 
	qm2_scf.F90 
	qm2_dftb_ewald.F90 
	qm2_dftb_scf.F90 )

set(COMMON_SOURCES 	findmask.F90 file_io_dat.F90 qm2_allocate_e_repul.F90 qm2_calc_charges.F90 qm2_calc_dipole.F90 
	qm2_calc_rij_and_eqns.F90 qm2_dihed.F90 qm2_fock.F90 
	qm2_get_qm_forces.F90 qm2_get_qmmm_forces.F90 qm2_h1elec.F90 qm2_hcore_qmqm.F90 
	qm2_hcore_qmmm.F90 qm2_identify_peptide_links.F90 
	qm2_load_params_and_allocate.F90 qm2_repp.F90 qm2_rotate_qmqm.F90 
	qm2_core_core_repulsion.F90 qm2_core_core_repulsion_dxyz.F90 
	qm2_iterator_mod.F90 qm2_diagonalizer_module.F90 
	qm2_setup_orb_exp.F90 qm2_smallest_number.F90 qm2_fock_predict.F90 
	qm_gb.F90 qm_zero_charges.F90 qm_assign_atom_types.F90 
	qm_link_atoms.F90 qm2_print_charges.F90 
	qmmm_qmtheorymodule.F90 qm2_print_bondorders.F90 
	qm2_pm6_hof_module.F90 qm2_params_module.F90
	qm2_dftb_module.F90 qm2_dftb_broyden.F90 
	qm2_dftb_dispersion_egr.F90 qm2_dftb_dispersion_params.F90 
	qm2_dftb_dispersionread.F90 qm2_dftb_energy.F90 
	qm2_dftb_ewevge.F90 qm2_dftb_externalshift.F90 
	qm2_dftb_fermi.F90 qm2_dftb_get_qm_forces.F90 qm2_dftb_gamma.F90 
	qm2_dftb_gammamat.F90 qm2_dftb_get_qmmm_forces.F90 qm2_dftb_gettab.F90 
	qm2_dftb_load_params.F90 
	qm2_dftb_mulliken.F90 
	qm2_dftb_repulsiv.F90 qm2_dftb_self.F90 qm2_dftb_shift.F90 
	qm2_dftb_skpar.F90 qm2_dftb_slkode.F90 
	qm2_dftb_slktrafo.F90 qm2_dftb_read_cm3.F90 qm2_dftb_cm3.F90 
	qm2_dftb_gb.F90 qm2_dftb_3rd_order.F90 
	qmmm_module.F90 qmmm_nml_module.F90 qmmm_struct_module.F90 qmmm_vsolv_module.F90 
	qm2_parameters.F90 dh_correction_module.F90 
	utilitiesModule.F90 elementOrbitalIndex.F90 parameterReader.F90 
	slater_overlap.F90 qm2_h1elec_d.F90 rotation.F90 qm2_repp_d.F90 
	opnq_Edisp.F90  opnq_Erep.F90  opnq_Evdw.F90  opnq.F90  opnq_SwitchMod.F90 
	qm2_fock_d.F90 MNDOChargeSeparation.F90 qm2_print_energy.F90
	qm2_dftb_get_dftb3_parameters.F90 qm2_dftb_gamma_dftb3.F90)

# --------------------------------------------------------------------

include_directories(${AMBER_COMMON_MOD_DIR})

#handle optimization flags
set(ALL_FORTRAN_SOURCES ${EXECUTABLE_SOURCES} ${LIBRARY_SOURCES} ${COMMON_SOURCES})
apply_optimization_declarations(${ALL_FORTRAN_SOURCES})

# QUICK support
if(BUILD_QUICK)
    add_definitions(-DQUICK)
endif()

# TCPB support
if(BUILD_TCPB)
    add_definitions(-DTCPB)
endif()

# --------------------------------------------------------------------

add_library(sqm_common OBJECT ${COMMON_SOURCES})
set_property(TARGET sqm_common PROPERTY POSITION_INDEPENDENT_CODE ${SHARED})
config_module_dirs(sqm_common ${SQMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
add_dependencies(sqm_common amber_common)

add_executable(sqm ${EXECUTABLE_SOURCES} $<TARGET_OBJECTS:sqm_common>)
target_link_libraries(sqm sff_fortran netlib amber_common)
set_property(TARGET sqm PROPERTY COMPILE_DEFINITIONS SQM)
config_module_dirs(sqm ${SQMEXE_MOD_DIR} ${SQMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

#NOTE: currently must be static, because it contains undefined references defined in the SANDER code
add_library(libsqm STATIC ${LIBRARY_SOURCES} $<TARGET_OBJECTS:sqm_common>)
config_module_dirs(libsqm ${SQMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
target_link_libraries(libsqm netlib amber_common)
make_pic_if_needed(libsqm)
install_libraries(libsqm)
remove_prefix(libsqm)


install(TARGETS sqm DESTINATION ${BINDIR})

# --------------------------------------------------------------------
# Parallel versions
# --------------------------------------------------------------------
if(MPI)
	make_mpi_version(sqm_common sqm_common_mpi LANGUAGES Fortran)
	config_module_dirs(sqm_common_mpi ${SQMLIB_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
	add_dependencies(sqm_common_mpi amber_common)
	
	make_mpi_version(sqm sqm.MPI LANGUAGES Fortran SWAP_SOURCES $<TARGET_OBJECTS:sqm_common> TO $<TARGET_OBJECTS:sqm_common_mpi> INSTALL)
	config_module_dirs(sqm.MPI ${SQMEXE_MPI_MOD_DIR} ${SQMLIB_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
	
	make_mpi_version(libsqm libsqm_mpi LANGUAGES Fortran SWAP_SOURCES $<TARGET_OBJECTS:sqm_common> TO $<TARGET_OBJECTS:sqm_common_mpi>)
	config_module_dirs(libsqm_mpi ${SQMLIB_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
endif()

if(OPENMP)
	make_openmp_version(sqm_common sqm_common_omp LANGUAGES Fortran)
	config_module_dirs(sqm_common_omp ${SQMLIB_OMP_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
	target_compile_definitions(sqm_common_omp PRIVATE OPENMP)
	add_dependencies(sqm_common_omp amber_common)

	make_openmp_version(libsqm libsqm_omp LANGUAGES Fortran SWAP_SOURCES $<TARGET_OBJECTS:sqm_common> TO $<TARGET_OBJECTS:sqm_common_omp>)
	config_module_dirs(libsqm_omp ${SQMLIB_OMP_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
	target_compile_definitions(libsqm_omp PRIVATE OPENMP)
endif()


if(xtb_EXTERNAL)
  add_compile_definitions(XTB)
endif()

if(dftbplus_EXTERNAL)
  add_compile_definitions(DFTBPLUS)
endif()
