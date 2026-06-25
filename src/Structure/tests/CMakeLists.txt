# Tests for qcStructure

add_executable(UTStructure ../../Common/tests/gtestmain.C
    StructureTests.C
    MolecularMesh1Tests.C
)

target_compile_options(UTStructure PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTStructure PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTStructure PROPERTY CXX_STANDARD 20)

target_include_directories(UTStructure PRIVATE
    # ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTStructure PRIVATE qcStructure qcMesh1 gtest)

gtest_discover_tests(UTStructure)
