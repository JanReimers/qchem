# Integrals-only microbenchmark (M&D vs libcint).  Not a unit test -- a profiling/perf driver.
add_executable(MnDBench MnDBench.C)
target_compile_definitions(MnDBench PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(MnDBench PRIVATE ${BasisDataPath})
set_property(TARGET MnDBench PROPERTY CXX_STANDARD 20)
target_include_directories(MnDBench PRIVATE
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_JSON_INCLUDE})
target_link_libraries(MnDBench PRIVATE qcGaussian_BS lapack)
