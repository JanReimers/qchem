# Tests for qcGaussian_BS

add_executable(UTGaussian_BS ../../Radial/tests/gtestmain.C
   M_PGSymmetry.C
   M_PG_BoxWalk.C
   M_LatticeScreener.C
   M_PG_MeshIntegrals.C
   M_PG_Oracle.C
   M_Evaluator.C
   M_MEvaluator.C
   M_Spherical.C
   M_LibCint.C
)

target_compile_definitions(UTGaussian_BS PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(UTGaussian_BS PRIVATE ${BasisDataPath})
target_compile_definitions(UTGaussian_BS PRIVATE ORACLE_DATA_PATH="${PROJECT_SOURCE_DIR}/tools/oracle")
set_property(TARGET UTGaussian_BS PROPERTY CXX_STANDARD 20)

target_include_directories(UTGaussian_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTGaussian_BS PRIVATE qcGaussian_BS lapack gtest)

gtest_discover_tests(UTGaussian_BS)
