add_library(qcLattice_BS) 
target_sources(qcLattice_BS PRIVATE Imp/Factory.C
PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES Factory.C
)
target_include_directories(qcLattice_BS PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
set_property(TARGET qcLattice_BS PROPERTY CXX_STANDARD 20)
target_link_libraries(qcLattice_BS PUBLIC qcBasisSet qcSymmetry qcCluster qcMesh qcCommon)
add_subdirectory(PolarizedGaussian)