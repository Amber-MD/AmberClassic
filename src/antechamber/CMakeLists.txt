project(antechamber C CXX)

if("${CMAKE_C_COMPILER_ID}" STREQUAL GNU)
	add_compile_options(-Wno-sign-compare)
endif()

add_compile_options(${OPT_CFLAGS})

# common eprintf error-handling object
add_library(antechamber_eprintf_obj OBJECT eprintf.c)

#most antechamber programs are compiled from a single .c file
# --------------------------------------------------------------------
set(PROGRAMS am1bcc atomtype espgen parmcal parmchk2 prepgen residuegen respgen match match_atomname)

foreach(PROGRAM ${PROGRAMS})
	add_executable(${PROGRAM} "${PROGRAM}.c" $<TARGET_OBJECTS:antechamber_eprintf_obj>)
	target_link_libraries(${PROGRAM} m) 
endforeach()

target_link_libraries(cifparse)

#now get the irregular ones
# --------------------------------------------------------------------
add_executable(antechamber antechamber.c $<TARGET_OBJECTS:antechamber_eprintf_obj>)
target_link_libraries(antechamber cifparse)

add_executable(bondtype bondtype.c ncsu-penalties.c $<TARGET_OBJECTS:antechamber_eprintf_obj>)

# add_executable(bondtype.orig bondtype.cpp $<TARGET_OBJECTS:antechamber_eprintf_obj>)

add_executable(test-penalties test-penalties.c ncsu-penalties.c $<TARGET_OBJECTS:antechamber_eprintf_obj>)

add_executable(test-int-part ncsu-penalties.c $<TARGET_OBJECTS:antechamber_eprintf_obj>)
set_target_properties(test-int-part PROPERTIES COMPILE_DEFINITIONS ENABLE_INT_PART_MAIN=1)

add_executable(test-combination ncsu-penalties.c $<TARGET_OBJECTS:antechamber_eprintf_obj>)
set_target_properties(test-combination PROPERTIES COMPILE_DEFINITIONS ENABLE_COMBINATION_MAIN=1)

# --------------------------------------------------------------------

set(ALL_PROGRAMS ${PROGRAMS} antechamber bondtype test-penalties test-int-part test-combination)

targets_link_libraries(${ALL_PROGRAMS} LIBRARIES C::Math)

install_executables_wrapped(TARGETS ${PROGRAMS} antechamber bondtype)
