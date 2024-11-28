set(MAKEDIST_RST_SOURCES makeDIST_RST.c)

set(MAKEANG_RST_SOURCES makeANG_RST.c)

set(NEF_TO_RST_SOURCES nef_to_RST.c)

set(DATA_FILES_TO_INSTALL map.DG-AMBER map.NEF-AMBER tordef.lib)

#I noticed that this directory has a bunch of scripts, none of which are currently being installed
#if you ever want scripts to be installed, put them in this variable.
set(SCRIPTS_TO_INSTALL makeCHIR_RST makeCSA_RST.na makeRIGID_RST makeDIP_RST.dna makeDIP_RST.protein)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#not sure whether this should be optimized, so I'm guessing yes
add_compile_options(${OPT_CFLAGS})

add_executable(makeDIST_RST ${MAKEDIST_RST_SOURCES})
target_link_libraries(makeDIST_RST m)

add_executable(makeANG_RST ${MAKEANG_RST_SOURCES})
target_link_libraries(makeANG_RST m)

add_executable(nef_to_RST ${NEF_TO_RST_SOURCES})
target_link_libraries(nef_to_RST cifparse m)

install(TARGETS makeDIST_RST makeANG_RST nef_to_RST DESTINATION ${BINDIR})
install(FILES ${DATA_FILES_TO_INSTALL} DESTINATION ${DATADIR} COMPONENT Data)
install(PROGRAMS ${SCRIPTS_TO_INSTALL} DESTINATION ${BINDIR})

