# Tests for qcStructure

add_executable(UTStructure ../../Common/tests/gtestmain.C
    StructureTests.C
    MolecularMeshTests.C
    SymmetrizeMeshUT.C
    BeckeMeshUT.C
    BravaisUT.C
    SiteGroupsUT.C      # SiteRotations vs SiteEnvironmentRotations on MnO AFM-II and diamond (+U increment 2, 2026-09-21)
    StructureDataUT.C   # qchem.StructureData (D-STRUCTDATA)
    CrystalFilesUT.C    # qchem.CrystalFiles: the import/export translators (round trips, a COD CIF, the refusals)
    ImpositionUT.C      # ResolveImposition: the run's symmetry policy (V-BS-SYMMETRY S1)
    PointSamplingUT.C   # the point set chooses how a Fourier series is sampled on it (the rho-tilde sampler)
    CacheKeysUT.C       # a cache key is exact, never a printed double (KeyString, Structure::ID; 2026-11-12)
)

target_compile_definitions(UTStructure PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS} STRUCTURE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
set_property(TARGET UTStructure PROPERTY CXX_STANDARD 20)

target_include_directories(UTStructure PRIVATE
    # ${QCHEM_BLAZE_DIR}
    ${QCHEM_JSON_INCLUDE}
    ${QCHEM_GTEST_INCLUDE}
    # ${QCHEM_TABULATE_INCLUDE}
)

target_link_libraries(UTStructure PRIVATE qcStructure qcMesh gtest)

qchem_discover_tests(UTStructure)
