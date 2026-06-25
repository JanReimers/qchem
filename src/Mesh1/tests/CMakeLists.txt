# Tests for qcMesh1

add_executable(UTMesh1 ../../Common/tests/gtestmain.C
    Mesh1Primitives.C
)

target_compile_options(UTMesh1 PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTMesh1 PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTMesh1 PROPERTY CXX_STANDARD 20)

target_include_directories(UTMesh1 PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)

target_link_libraries(UTMesh1 PRIVATE qcMesh1 gtest)

gtest_discover_tests(UTMesh1)
