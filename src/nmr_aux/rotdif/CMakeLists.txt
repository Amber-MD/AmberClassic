
add_compile_options(${OPT_FFLAGS})

add_executable(randvec randvec.F90)

add_executable(rmscorr rmscorr.F90)

add_executable(tensorfit tensorfit.F90)

install(TARGETS randvec rmscorr tensorfit DESTINATION ${BINDIR})