# Tests for qcMesh

add_executable(UTMesh ../../Common/tests/gtestmain.C
    MeshPrimitives.C
)

target_compile_options(UTMesh PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTMesh PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTMesh PROPERTY CXX_STANDARD 20)

target_include_directories(UTMesh PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)

target_link_libraries(UTMesh PRIVATE qcMesh gtest)

gtest_discover_tests(UTMesh)
