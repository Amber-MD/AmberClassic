#CMake buildfile for Small Angle X-ray Scattering
project(saxs CXX)

#Source files common to all saxs executables and libraries 
set(SAXS_COMMON_SOURCES lebedev.cpp atm_f.cpp dx.cpp f_calc.cpp IO.cpp pdb.cpp)

#Source files in the saxs_rism executable and library
set(SAXS_RISM_SOURCES saxs_rism.cpp vol_decomp.cpp)

#source files in the saxs_md executable and library
set(SAXS_MD_SOURCES saxs_md.cpp)

#------------------------------------------------------------------------------------------

add_compile_options(${OPT_CXXFLAGS})

#add the common object libraries
add_library(saxs_common_obj OBJECT ${SAXS_COMMON_SOURCES})
make_pic_if_needed(saxs_common_obj)

#executables
add_executable(saxs_rism ${SAXS_RISM_SOURCES} $<TARGET_OBJECTS:saxs_common_obj>)
add_executable(saxs_md ${SAXS_MD_SOURCES} $<TARGET_OBJECTS:saxs_common_obj>)

install(TARGETS saxs_rism saxs_md DESTINATION ${BINDIR})

#libraries
add_library(libsaxs_rism  ${SAXS_RISM_SOURCES} $<TARGET_OBJECTS:saxs_common_obj>)
add_library(libsaxs_md  ${SAXS_MD_SOURCES} $<TARGET_OBJECTS:saxs_common_obj>)
remove_prefix(libsaxs_rism)
remove_prefix(libsaxs_md)

install_libraries(libsaxs_rism libsaxs_md)



if(OPENMP)
	
	# check for an OpenMP implementation that supports addition in a #pragma omp atomic
	
	set(CMAKE_REQUIRED_FLAGS ${OpenMP_C_FLAGS})
	
	check_c_source_compiles(

"int main()
{
	int foo = 2;	
	#pragma omp atomic
		foo = foo * 2;
	return 0;
}"
		 OMP_ATOMIC_ADDITION_SUPPORT)
	
	if(NOT OMP_ATOMIC_ADDITION_SUPPORT)
		message(WARNING "Your compiler's OpenMP implementation is too old to compile saxs. OpenMP version 3.1 (spec date 201307) is required.  saxs will be built as serial only.")
	else()
		make_openmp_version(saxs_common_obj saxs_common_obj_omp LANGUAGES C)
		
		make_openmp_version(saxs_rism saxs_rism.OMP LANGUAGES C SWAP_SOURCES $<TARGET_OBJECTS:saxs_common_obj> TO $<TARGET_OBJECTS:saxs_common_obj_omp> INSTALL)
		make_openmp_version(saxs_md saxs_md.OMP LANGUAGES C SWAP_SOURCES $<TARGET_OBJECTS:saxs_common_obj> TO $<TARGET_OBJECTS:saxs_common_obj_omp> INSTALL)
		
		make_openmp_version(libsaxs_rism libsaxs_rism.OMP LANGUAGES C SWAP_SOURCES $<TARGET_OBJECTS:saxs_common_obj> TO $<TARGET_OBJECTS:saxs_common_obj_omp> INSTALL)
		make_openmp_version(libsaxs_md libsaxs_md.OMP LANGUAGES C SWAP_SOURCES $<TARGET_OBJECTS:saxs_common_obj> TO $<TARGET_OBJECTS:saxs_common_obj_omp> INSTALL)
	endif()
endif()
