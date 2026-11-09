# Tests for qcGaussian_BS

add_executable(UTGaussian_BS ../../Radial/tests/gtestmain.C
   M_PGSymmetry.C
   M_SALCRepresentation.C   # a molecular SALC block end to end: Identity exact against the raw block's Dense rep (V-BS-SYMMETRY S3)
   M_PG_BoxWalk.C
   M_TransitionCollocation.C   # the (k+q, k) transition collocation + adjoint (LinearResponsePlan B1, R3 step 2)
   M_LatticeScreener.C
   M_PG_MeshIntegrals.C
   M_Gaussian94Reader.C   # a missing element THROWS naming element + file (was a garbage-byte parse, 2026-10-21)
   M_PG_Oracle.C
   M_Evaluator.C
   M_MEvaluator.C
   M_Spherical.C
   M_GridResolution.C   # the grid builders' question (GridResolution + tBasisSet::MakeResolution; V-BS-RASTER, 2026-11-02)
   M_PP_ClosedForm.C    # the closed-form molecular PP integrals against a fine mesh
   M_BareCoulomb.C      # the BareCoulombSource face vs the d-shell Slater averages (ACBN0, +U increment 3; 2026-09-21)
   M_LibCint.C
   M_GeometryDerivatives.C   # the derivative methods on the integral faces: PointJets, the PG_Cart pullbacks, the refusals (§13.2)
   DenseIntegratorUT.C
)

target_compile_definitions(UTGaussian_BS PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(UTGaussian_BS PRIVATE ${BasisDataPath})
target_compile_definitions(UTGaussian_BS PRIVATE ORACLE_DATA_PATH="${PROJECT_SOURCE_DIR}/tools/oracle")
set_property(TARGET UTGaussian_BS PROPERTY CXX_STANDARD 20)

target_include_directories(UTGaussian_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/src   # forward.H (PointJetsTests)
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTGaussian_BS PRIVATE qcGaussian_BS lapack gtest)

qchem_discover_tests(UTGaussian_BS)
