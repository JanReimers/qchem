# Tests for qcResponse (doc/LinearResponsePlan.md).  Model-driven: the Reference and the probes are built
# from VALUES, so the sum-over-states machinery is checked against an exactly solvable brute-force model
# with no SCF in the loop.

add_executable(UTResponse ../../Common/tests/gtestmain.C
    Ring.C              # chi0(q) -> chi0(R) vs a brute-force finite ring (insulator + metal)
    Kernel.C            # the response kernel vs its finite-difference oracle (D6), on H2O
)

target_compile_definitions(UTResponse PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTResponse PROPERTY CXX_STANDARD 20)

target_include_directories(UTResponse PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/src                       # forward.H (the FD oracle's friend)
)

target_link_libraries(UTResponse PRIVATE qcResponse qchem lapack gtest)   # qchem: the Gaussian basis factory (the H2O kernel gates)

qchem_discover_tests(UTResponse)
