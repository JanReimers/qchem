# Tests for qcMolecule_BS

add_executable(UTMolecule_BS ../../Atom/tests/gtestmain.C
   M_PGSymmetry.C 
   M_PG_MeshIntegrals.C
)

target_compile_options(UTMolecule_BS PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTMolecule_BS PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS} BASISSET_DATA_PATH=\"${PROJECT_SOURCE_DIR}/BasisSetData\")
set_property(TARGET UTMolecule_BS PROPERTY CXX_STANDARD 20)

target_include_directories(UTMolecule_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/wignerSymbols/include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTMolecule_BS PRIVATE qcMolecule_BS lapack gtest)

gtest_discover_tests(UTMolecule_BS)
