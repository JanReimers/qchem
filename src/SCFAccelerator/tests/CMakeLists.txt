# Tests for qcSCFAccelerator

add_executable(UTSCFAccelerator ../../Common/tests/gtestmain.C
    SCFAccelerator.C
)

target_compile_definitions(UTSCFAccelerator PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTSCFAccelerator PROPERTY CXX_STANDARD 20)

target_include_directories(UTSCFAccelerator PRIVATE
    ${QCHEM_GTEST_INCLUDE}
    ${QCHEM_BLAZE_DIR}
)

target_link_libraries(UTSCFAccelerator PRIVATE qcSCFAccelerator lapack gtest)

qchem_discover_tests(UTSCFAccelerator)
