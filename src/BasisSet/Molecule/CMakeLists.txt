add_library(qcMolecule_BS)
target_sources(qcMolecule_BS PRIVATE Imp/Factory.C
Imp/SymmetryAdaptedBasisSet.C
Imp/BasisFiles.C
Readers/Imp/Gaussian94.C
PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES Factory.C
SymmetryAdaptedBasisSet.C
IrrepBasisSet.C
LatticeSum1E.C
BasisFiles.C
Reader.C
Readers/Gaussian94.C
)
target_include_directories(qcMolecule_BS PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
set_property(TARGET qcMolecule_BS PROPERTY CXX_STANDARD 20)
target_compile_options(qcMolecule_BS PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcMolecule_BS PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(qcMolecule_BS PRIVATE ${BasisDataPath})

target_link_libraries(qcMolecule_BS PUBLIC qcBasisSet qcElConfig qcStructure)
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
add_subdirectory(PG_Cart)
add_subdirectory(PG_Spherical)
add_subdirectory(PG_LibCint)
add_subdirectory(Evaluators)
add_subdirectory(tests)
add_subdirectory(bench)