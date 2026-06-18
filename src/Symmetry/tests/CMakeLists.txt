# Tests for qcSymmetry

add_executable(UTSymmetry ../../Common/tests/gtestmain.C
   SymQN.C
   M_PointGroup.C 
   M_CartesianRep.C 
   M_CharacterTable.C 
   M_SALC.C
)

target_compile_options(UTSymmetry PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTSymmetry PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTSymmetry PROPERTY CXX_STANDARD 20)

target_include_directories(UTSymmetry PRIVATE
    # ${PROJECT_SOURCE_DIR}/submodules/blaze
    # ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTSymmetry PRIVATE qcSymmetry gtest)

gtest_discover_tests(UTSymmetry)
