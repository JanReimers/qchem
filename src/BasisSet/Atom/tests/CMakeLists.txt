# Tests for qcAtom_BS

add_executable(UTAtom_BS gtestmain.C
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

target_compile_options(UTAtom_BS PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTAtom_BS PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTAtom_BS PROPERTY CXX_STANDARD 20)

target_include_directories(UTAtom_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTAtom_BS PRIVATE qcAtom_BS qcLASolver gtest)

gtest_discover_tests(UTAtom_BS)
