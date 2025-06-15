
file (GLOB_RECURSE qchem_SOURCES CONFIGURE_DEPENDS "*.C")
file (GLOB_RECURSE qchem_FORTRAN CONFIGURE_DEPENDS "*.f")
list(FILTER qchem_SOURCES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/Common/.*" )
list(FILTER qchem_SOURCES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/Symmetry/.*" )
#message("qchem_SOURCES is ${qchem_SOURCES}")
target_sources(qchem PRIVATE ${qchem_SOURCES} ${qchem_FORTRAN})

add_subdirectory(Common)
add_subdirectory(Symmetry)
