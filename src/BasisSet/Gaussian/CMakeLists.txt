# src/BasisSet/Gaussian — the GAUSSIAN engine library (doc/BasisSetTaxonomyPlan.md §3; was Molecule/ until
# step 2, 2026-09-13).  FAMILY = Gaussian {primitive, contracted} x {Cartesian, spherical}; the G is the
# MODULE tag, not the library:  Gaussian.Evaluators.*  the engine (no G -- serves both);
#                               Gaussian.Point.*       G = P, the molecular bases (Point/);
#                               Gaussian.Lattice.*     G = T, the Bloch-sum half (Lattice/: LatticeSum1E,
#                                                      LatticeScreener, the spherical lattice view, GPW).
# Invariant (plan step 3, ctest audit): no .Point. module imports qchem.UnitCell / qchem.Lattice_3D /
# qchem.Symmetry.Lattice_3D.*.
add_library(qcGaussian_BS)
# The module file set is declared HERE (base dir = this directory) so the Evaluators/ Point/ Lattice/
# subdirectories can all contribute to it; the first declaration fixes the base directory.
target_sources(qcGaussian_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(qcGaussian_BS PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
set_property(TARGET qcGaussian_BS PROPERTY CXX_STANDARD 20)
target_compile_definitions(qcGaussian_BS PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(qcGaussian_BS PRIVATE ${BasisDataPath})

# qcPlaneWave_BS: the GPW engine is (Gaussian orbitals) x (PW density fit) -- it holds a PW_Grid_Evaluator
# and returns PlaneWaveFit_IBS fit bases.  A frozen family pairing (plan §5); the edge runs THIS way only.
target_link_libraries(qcGaussian_BS PUBLIC qcBasisSet qcElConfig qcStructure qcPlaneWave_BS)
# OpenMP for the GPW pair loops (PG_Cart_MnD evaluator) is applied GLOBALLY at the top level -- C++20
# modules require one OpenMP language mode across the whole import graph (a non-OpenMP PCM cannot be
# imported by an -fopenmp TU), so the old per-target scoping here cannot build; see the top-level
# CMakeLists OpenMP block for the full story (incl. the 2026-07-31 finding that clang's
# -fopenmp=libgomp fallback generated NO parallel code and ran silently serial).
if(NOT OpenMP_CXX_FOUND)
    message(WARNING "OpenMP (libomp) not found: GPW_OMP_THREADS will be SILENTLY SERIAL.  clang cannot "
                    "generate OpenMP code against libgomp.  Install libomp (sudo apt install libomp-dev) "
                    "and point cmake at it (see the top-level CMakeLists OpenMP block).")
endif()
add_subdirectory(Evaluators)
add_subdirectory(Point)
add_subdirectory(Lattice)
add_subdirectory(tests)
add_subdirectory(bench)