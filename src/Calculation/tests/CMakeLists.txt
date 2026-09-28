# Tests for the qchem::Calculation facades

add_executable(UTCalculation ../../Common/tests/gtestmain.C
    RunDiagnostics.C
    Materials.C
    ValenceBasisGen.C   # the valence-basis generator (from IntegrationTests 2026-09-15, TE phase 5)
    SolidState.C        # saved periodic SCF states: save / exact resume / warm start / refusal (CK-1)
)

target_compile_definitions(UTCalculation PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTCalculation PROPERTY CXX_STANDARD 20)

target_include_directories(UTCalculation PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)

target_link_libraries(UTCalculation PRIVATE qchem gtest)

qchem_discover_tests(UTCalculation)
