# Tests for qcSymmetry

add_executable(UTSymmetry ../../Common/tests/gtestmain.C
   SymQN.C
   M_PointGroup.C 
   M_CartesianRep.C
   M_SphericalRep.C
   M_CharacterTable.C
   M_SALC.C
   L_SpaceGroup.C
   L_Fold.C
   RepresentationUT.C   # Representation<T> + the spin-pair projector on hand-built reps (V-BS-SYMMETRY S3)
)

target_compile_definitions(UTSymmetry PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTSymmetry PROPERTY CXX_STANDARD 20)

target_include_directories(UTSymmetry PRIVATE
    # ${QCHEM_BLAZE_DIR}
    # ${QCHEM_JSON_INCLUDE}
    # ${QCHEM_TABULATE_INCLUDE}
    ${QCHEM_GTEST_INCLUDE}
)

target_link_libraries(UTSymmetry PRIVATE qcSymmetry gtest)

qchem_discover_tests(UTSymmetry)
