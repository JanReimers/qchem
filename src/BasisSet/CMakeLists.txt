add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
  Imp/BasisSet.C
  Imp/Conditioning.C
  Imp/BareCoulombSource.C
  Imp/AtomCentredSet.C
  Imp/PointImages.C
  Internal/Imp/Fit_IBS.C
  Imp/DeltaFit_IBS.C
  Internal/Imp/SymmetryAdapted_IBS.C
  Internal/Imp/Orbital_DHF_IBS.C
  Internal/Imp/Orbital_ERI4_IBS.C
  Internal/Imp/CacheSlot.C
  Internal/Imp/CachedIntegrals.C
  Internal/Imp/Cache4.C
  Internal/Imp/Cache2.C
  Internal/Imp/Cache3.C
  Imp/PairProjection.C
  Imp/GSpaceDensityProjector.C
  Imp/RasterDensityProjector.C
  Internal/Imp/ERI4.C
  Internal/Imp/ERI4T.C

)
target_sources(qcBasisSet PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    DBCacheClient.C                       # the integral cache's KEY CONTRACT (V1.20d; the mechanism stays Internal)
    Orbital_PP_IBS.C                      # the species-field INTEGRAL SERVICE + its ARGUMENT vocabulary (V1.2; absorbed SpeciesField.C 2026-09-15)
    PointSetPullback.C                    # a block's point-set quadrature pullback (the δ XC force on the Becke mesh; lane B 2026-11-19)
  BasisSet.C
  Conditioning.C                        # the basis conditioning report, one phase on every facade (pin 7)
  GridResolution.C                      # 'how fine and how far must a grid over me reach' -- the grid builders' query (V-BS-RASTER; 2026-11-02)
  AoShellSource.C                       # 'I am built from atom-centred shells' (+U manifolds, SALC induction; 2026-09-19)
  PointImages.C                         # the k-independent lattice images of a Bloch basis on a point set (shared Φ; 2026-11-19)
  OrbitalImageSource.C                  # 'the image of my Bloch orbitals under a crystal op' (KAxisPlan P4, D5/CK-3a; 2026-10-27)
  AtomCentredSet.C                      # an atom-centred function set + 'my overlap with it' (+U projectors, HubbardProjectorPlan D1/D2)
  BareCoulombSource.C                   # 'bare (m1m2|m3m4) over a function subset' -- the ACBN0 face (+U increment 3; 2026-09-21)
  IrrepBasisSet.C
  PointJets.C                           # :PointJets, a partition of qchem.BasisSet.IrrepBasisSet: the slab of function jets (§13.2; 2026-11-08)
  Orbital_1E_IBS.C
  ImplicitAngular_IBS.C
  G_FieldEvaluator.C
  G_RasterStarField.C                   # a raster field the density group star-averages (V-SYM-SEQUESTER-OPS; 2026-11-09)
  Orbital_HF_IBS.C
  Orbital_DFT_IBS.C
  Transition_DFT_IBS.C                  # the (k+q, k) transition-density faces (LinearResponsePlan B2)
  Orbital_DHF_IBS.C
    Lattice_IBS.C                         # the G=T SPEC tier: isLattice_{1E,DFT}_Evaluator + the Lattice::*_IBS<E,T> mixins (TaxonomyPlan 1a0)
  Fit_Types.C
  Fit_IBS.C
  Internal/Fit_IBS.C
  DeltaFit_IBS.C
  GradOverlap3C.C                       # a point fit basis's gradient 3-centre overlap (GGA stage 2; 2026-11-21)
  MultipoleFit_IBS.C
  CoulombDeltaFit_IBS.C
  GSpaceDensityProjector.C              # the ρ̃ projector FACTORY (returns the PairProjector<dcmplx> face)
  Internal/GSpaceDensityProjector.C     # ...its class (only the factory names it)
  RasterDensityProjector.C              # the raw rho_DM projector on a fit raster (R1.0n/o, 2026-11-19)
  Internal/SymmetryAdapted_IBS.C
  Internal/BasisSetImp.C
  Internal/IrrepBasisSetImp.C
  Internal/CacheSlot.C                  # a block's slot in its parent's cache + the cached overlap (D-BASISCACHE-FANIN, 2026-11-18)
  Internal/CachedIntegrals.C            # the cached kinetic and 3-centre accessors through it
  BasisCache.C                          # the parent-made cache of a basis's own tables (IntegralCacheOwnershipPlan S3a)
  Internal/Orbital_DHF_IBS.C
  Internal/Orbital_ERI4_IBS.C
  Internal/Cache4.C
  Internal/Cache2.C
  Internal/Cache3.C
  Internal/PrimitiveStore.C             # one owner's Ω / RNLM / H3 caches (IntegralCacheOwnershipPlan S5)
  PairProjection.C
  PairProjector.C
  GMap.C
  Internal/ERI4.C
  Internal/ERI4T.C
  Internal/IntegralEnums.C
  )
target_include_directories(qcBasisSet PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_property(TARGET qcBasisSet PROPERTY CXX_STANDARD 20)
target_include_directories(qcBasisSet PRIVATE ${QCHEM_BLAZE_DIR})
target_link_libraries(qcBasisSet PUBLIC qcElConfig qcStructure qcMesh qcLASolver qcAD)   # qcAD: NoGradient, the derivative faces' refusal
target_compile_definitions(qcBasisSet PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
add_subdirectory(Radial)      # the O(3) engine (was Atom/; TaxonomyPlan step 2)
add_subdirectory(PlaneWave)   # the plane-wave engine: a leaf over the core (TaxonomyPlan 1a)
add_subdirectory(Gaussian)    # the Gaussian engine: Point/ (G=P) + Lattice/ (G=T) tiers (TaxonomyPlan step 2)
add_subdirectory(Lattice)     # the T CONTAINER tier over both engines (was Lattice_3D/; TaxonomyPlan step 2)

