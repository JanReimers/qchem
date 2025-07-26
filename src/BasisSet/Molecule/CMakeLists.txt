add_library(qcMolecule_BS) 
target_sources(qcMolecule_BS PRIVATE Imp/Factory.C
PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES Factory.C
)
target_include_directories(qcMolecule_BS PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
set_property(TARGET qcMolecule_BS PROPERTY CXX_STANDARD 20)
target_link_libraries(qcMolecule_BS PUBLIC qcBasisSet qcSymmetry qcCluster qcMesh qcCommon)
add_subdirectory(PolarizedGaussian)