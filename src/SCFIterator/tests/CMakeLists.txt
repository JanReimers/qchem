# Tests for the SCFIterator.  These drive a real (small) SCF through AtomCalculation -- the iterator
# cannot be constructed without a basis/Hamiltonian/accelerator -- but what they ASSERT is the iterator's
# own behaviour (the trace column set), so they belong here rather than in IntegrationTests.
# Linking `qchem` is not a layering inversion: SCFIterator and Calculation are both part of that ONE
# aggregate library, so this exe pulls in no library that sits above the code under test.

add_executable(UTSCFIterator ../../Common/tests/gtestmain.C
    SCFTrace.C
    MixingMode.C
    TaylorEnergy.C
)

target_compile_definitions(UTSCFIterator PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTSCFIterator PROPERTY CXX_STANDARD 20)

target_include_directories(UTSCFIterator PRIVATE
    ${QCHEM_GTEST_INCLUDE}
    ${QCHEM_BLAZE_DIR}
)

target_link_libraries(UTSCFIterator PRIVATE qchem lapack gtest)

qchem_discover_tests(UTSCFIterator)
