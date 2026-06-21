# Integrals-only microbenchmark (M&D vs libcint).  Not a unit test -- a profiling/perf driver.
add_executable(MnDBench MnDBench.C)
target_compile_options(MnDBench PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(MnDBench PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(MnDBench PRIVATE ${BasisDataPath})
set_property(TARGET MnDBench PROPERTY CXX_STANDARD 20)
target_include_directories(MnDBench PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include)
target_link_libraries(MnDBench PRIVATE qcMolecule_BS lapack)
