# Command-line applications (not tests).

# scfrun: command-line SCF driver for accelerator tuning (no recompile).  Drives qchem::AtomCalculation
# directly (no QchemTester); needs the TestUtils module for the RelativeHF/DHFError oracle helpers.
add_executable(scfrun scfrun.C)
target_sources(scfrun PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES
    BASE_DIRS ${PROJECT_SOURCE_DIR} FILES ../IntegrationTests/TestUtils.C)
target_compile_options(scfrun PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(scfrun PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET scfrun PROPERTY CXX_STANDARD 20)
target_include_directories(scfrun PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)
target_link_libraries(scfrun PUBLIC qchem qcMath qcCommon blas)
