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
# OpenMP for the GPW pair loops (PG_Cart_MnD evaluator).  Scoped to THIS lib's own TUs, and QCHEM_OPENMP
# is defined ONLY when a real OpenMP (libomp) is found -- so the code's parallel branch exists exactly
# when the pragmas actually codegen.
# MEASURED 2026-07-31: the old fallback `-fopenmp=libgomp` is GCC LINK-COMPATIBILITY ONLY in clang -- it
# links libgomp but generates NO parallel code (zero GOMP_*/omp_outlined symbols in every object; 99% CPU
# under GPW_OMP_THREADS=8; a hello-world OMP program compiled the same way runs serial).  So the fallback
# ran SILENTLY SERIAL since it was added.  clang needs libomp for codegen:
#     sudo apt install libomp-dev     (libomp.so + omp.h; the __kmpc ABI is version-stable)
# then re-run cmake -- find_package(OpenMP) succeeds and the branch below activates.
if(OpenMP_CXX_FOUND)
    target_link_libraries(qcMolecule_BS PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(qcMolecule_BS PRIVATE QCHEM_OPENMP)
else()
    message(WARNING "OpenMP (libomp) not found: GPW_OMP_THREADS will be SILENTLY SERIAL.  clang cannot "
                    "generate OpenMP code against libgomp.  Install libomp (sudo apt install libomp-dev) "
                    "and re-run cmake.")
endif()
add_subdirectory(PG_Cart)
add_subdirectory(PG_Spherical)
add_subdirectory(PG_LibCint)
add_subdirectory(Evaluators)
add_subdirectory(tests)
add_subdirectory(bench)