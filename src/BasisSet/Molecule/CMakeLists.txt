add_library(qcMolecule_BS)
target_sources(qcMolecule_BS PRIVATE Imp/Factory.C
Imp/SymmetryAdaptedBasisSet.C
Imp/BasisFiles.C
Readers/Imp/Gaussian94.C
PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES Factory.C
SymmetryAdaptedBasisSet.C
AoShellSource.C
IrrepBasisSet.C
BasisFiles.C
Reader.C
Readers/Gaussian94.C
)
target_include_directories(qcMolecule_BS PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
set_property(TARGET qcMolecule_BS PROPERTY CXX_STANDARD 20)
target_compile_options(qcMolecule_BS PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcMolecule_BS PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(qcMolecule_BS PRIVATE ${BasisDataPath})

target_link_libraries(qcMolecule_BS PUBLIC qcBasisSet qcElConfig qcStructure)
add_subdirectory(PG_Cart)
add_subdirectory(PG_Spherical)
add_subdirectory(PG_LibCint)
add_subdirectory(Evaluators)
add_subdirectory(tests)
add_subdirectory(bench)