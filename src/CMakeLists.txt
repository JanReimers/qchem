
file (GLOB_RECURSE qchem_SOURCES CONFIGURE_DEPENDS "*.C")
file (GLOB_RECURSE qchem_FORTRAN CONFIGURE_DEPENDS "*.f")
list(FILTER qchem_SOURCES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/Common/.*" )
list(FILTER qchem_SOURCES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/Symmetry/.*" )
list(FILTER qchem_SOURCES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/Mesh/.*" )
list(FILTER qchem_SOURCES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/Cluster/.*" )
list(FILTER qchem_SOURCES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/BasisSet/.*" )
list(FILTER qchem_SOURCES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/LASolver/.*" )
#message("qchem_SOURCES is ${qchem_SOURCES}")
target_sources(qchem PRIVATE ${qchem_SOURCES} ${qchem_FORTRAN})

add_subdirectory(Common)
add_subdirectory(Symmetry)
add_subdirectory(Mesh)
add_subdirectory(Cluster)
add_subdirectory(LASolver)
add_subdirectory(BasisSet)
