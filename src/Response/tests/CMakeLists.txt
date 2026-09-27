# Tests for qcResponse (doc/LinearResponsePlan.md).  Model-driven: the Reference and the probes are built
# from VALUES, so the sum-over-states machinery is checked against an exactly solvable brute-force model
# with no SCF in the loop.

add_executable(UTResponse ../../Common/tests/gtestmain.C
    Ring.C              # chi0(q) -> chi0(R) vs a brute-force finite ring (insulator + metal)
)

target_compile_definitions(UTResponse PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTResponse PROPERTY CXX_STANDARD 20)

target_include_directories(UTResponse PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)

target_link_libraries(UTResponse PRIVATE qcResponse lapack gtest)

qchem_discover_tests(UTResponse)
