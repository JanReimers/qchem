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
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googlemock/include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTRadial_BS PRIVATE qcRadial_BS qcLASolver gtest)

gtest_discover_tests(UTRadial_BS)
