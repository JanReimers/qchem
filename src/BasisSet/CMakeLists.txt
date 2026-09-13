add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
    Imp/IntegralsCacheReport.C
  Imp/BasisSet.C
  Imp/IrrepBasisSet.C
  Imp/Orbital_1E_IBS.C
  Imp/Orbital_DFT_IBS.C
  Imp/Orbital_DHF_IBS.C
  Internal/Imp/Fit_IBS.C
  Imp/DeltaFit_IBS.C
  Internal/Imp/SymmetryAdapted_IBS.C
  Internal/Imp/Orbital_DHF_IBS.C
  Internal/Imp/Orbital_ERI4_IBS.C
  Internal/Imp/DB_Cache.C
  Internal/Imp/Cache4.C
  Internal/Imp/Cache2.C
  Internal/Imp/Cache3.C
  Imp/Projector3.C
  Internal/Imp/ERI4.C
  Internal/Imp/ERI4T.C
  Internal/Imp/DB_Cache_RAM.C
  Internal/Imp/DB_Cache_RAM_Report.C

)
target_sources(qcBasisSet PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    DBCacheClient.C                       # the integral cache's KEY CONTRACT (V1.20d; the mechanism stays Internal)
    IntegralsCacheReport.C                # the one orchestrator-side ask of the cache (V1.20d)
    SpeciesField.C                        # the species-field ARGUMENT vocabulary (V1.2)
    Orbital_PP_IBS.C                      # the species-field INTEGRAL SERVICE (V1.2; was Pseudopotential::Integrals_Pseudo)
  BasisSet.C
  IrrepBasisSet.C
  Orbital_1E_IBS.C
  ImplicitAngular_IBS.C
  G_FieldEvaluator.C
  Orbital_HF_IBS.C
  Orbital_DFT_IBS.C
  Orbital_DHF_IBS.C
    Lattice_IBS.C                         # the G=T SPEC tier: isLattice_{1E,DFT}_Evaluator + the Lattice::*_IBS<E,T> mixins (TaxonomyPlan 1a0)
  Fit_Types.C
  Fit_IBS.C
  Internal/Fit_IBS.C
  DeltaFit_IBS.C
  Internal/SymmetryAdapted_IBS.C
  Internal/BasisSetImp.C
  Internal/IrrepBasisSetImp.C
  Internal/Orbital_DHF_IBS.C
  Internal/Orbital_ERI4_IBS.C
  Internal/Cache4.C
  Internal/Cache2.C
  Internal/Cache3.C
  Projector3.C
  GMap.C
  Internal/ERI4.C
  Internal/ERI4T.C
  Internal/IntegralEnums.C
  Internal/DB_Cache.C
  Internal/DB_Cache_RAM.C
  )
target_include_directories(qcBasisSet PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_property(TARGET qcBasisSet PROPERTY CXX_STANDARD 20)
target_include_directories(qcBasisSet PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_link_libraries(qcBasisSet PUBLIC qcElConfig qcStructure qcMesh)
target_compile_definitions(qcBasisSet PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
add_subdirectory(Radial)      # the O(3) engine (was Atom/; TaxonomyPlan step 2)
add_subdirectory(PlaneWave)   # the plane-wave engine: a leaf over the core (TaxonomyPlan 1a)
add_subdirectory(Molecule)
add_subdirectory(Lattice_3D)

