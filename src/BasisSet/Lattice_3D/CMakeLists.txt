# src/BasisSet/Lattice_3D — plane-wave (3D periodic) basis set.
# See doc/PlaneWavePlan.md.
add_library(qcLattice_BS)
target_sources(qcLattice_BS PRIVATE
    Imp/PlaneWave_IBS.C
    Imp/APW_IBS.C
    Imp/LAPW_IBS.C
    Imp/BasisSet.C
    Imp/GTH_Potentials.C
)
target_sources(qcLattice_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Internal/GVectors.C
    Internal/KPlusG.C
    GTH_Potentials.C
    PlaneWave_IBS.C
    APW_IBS.C
    LAPW_IBS.C
    BandStructure.C
    BasisSet.C
)
target_include_directories(qcLattice_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/src
)
set_property(TARGET qcLattice_BS PROPERTY CXX_STANDARD 20)
target_compile_options(qcLattice_BS PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcLattice_BS PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS} ${LatticeDataPath})
target_link_libraries(qcLattice_BS PUBLIC qcBasisSet qcStructure qcSymmetry qcLASolver)

add_subdirectory(tests)
