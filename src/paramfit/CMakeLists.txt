project(paramfit C)

set(PARAMFIT_SOURCE calc_r_squared.c create_input.c defaults.c 
	elements.c error_messages.c eval_amber_std.c 
	file_io.c fitting_control.c help_functions.c 
	mem_alloc.c misc_utils.c options_summary.c 
	parameter_optimiser.c param_summary.c print_program_info.c 
	process_command_line.c process_job_control_setting.c 
	process_prmtop.c read_prmtop.c simplex.c write_input.c 
	genetic_algorithm.c bounds_check.c eval_amber_forces.c 
	dihedral_fitting.c wizard.c read_mdcrd.c read_energy.c)

# --------------------------------------------------------------------

add_compile_options(${OPT_CFLAGS})

# check for backtrace function
check_all_symbols(execinfo.h backtrace)

if(NOT HAVE_BACKTRACE)
	add_definitions(-DWITHOUT_BACKTRACE)
endif()

# --------------------------------------------------------------------

add_executable(paramfit ${PARAMFIT_SOURCE})

#Add the math library if it exists

target_link_libraries(paramfit C::Math)

install(TARGETS paramfit DESTINATION ${BINDIR})

# --------------------------------------------------------------------
# Parallel Versions
# --------------------------------------------------------------------

if(OPENMP)
	make_openmp_version(paramfit paramfit.OMP LANGUAGES C INSTALL)
endif()