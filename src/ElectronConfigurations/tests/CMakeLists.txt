# Tests for qcElConfig

add_executable(UTElConfig ../../Common/tests/gtestmain.C
   ElectronConfigurations.C
   ElectronConfigurations_Dirac.C
)

target_compile_options(UTElConfig PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTElConfig PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTElConfig PROPERTY CXX_STANDARD 20)

target_include_directories(UTElConfig PRIVATE
    # ${PROJECT_SOURCE_DIR}/submodules/blaze
    # ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)

target_link_libraries(UTElConfig PRIVATE qcElConfig gtest)

gtest_discover_tests(UTElConfig)
