# Tests for qcRadial_BS

add_executable(UTRadial_BS gtestmain.C
   A_Cache4.C
   A_Pool.C
   Angular.C 
   BasisSet_Atom.C
   BasisSet_DHF.C 
   BSplines.C 
   DBCache.C
#    DiracHydrogen.C
   DiracIntegral.C 
   ERI4.C 
   GaussianIntegral.C 
#    GaussianERIIntegral.C
   RelAngular.C 
#    OpenMP.C
   SlaterIntegral.C 
)

target_compile_definitions(UTRadial_BS PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTRadial_BS PROPERTY CXX_STANDARD 20)

target_include_directories(UTRadial_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/src   # forward.H (A_Cache4's unit-test friend)
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_JSON_INCLUDE}
    ${QCHEM_BSPLINE_INCLUDE}
    ${QCHEM_GTEST_INCLUDE}
    ${googletest_SOURCE_DIR}/googlemock/include
    # ${QCHEM_TABULATE_INCLUDE}
)

target_link_libraries(UTRadial_BS PRIVATE qcRadial_BS qcLASolver gtest)

qchem_discover_tests(UTRadial_BS)
