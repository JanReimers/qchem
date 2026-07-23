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
# OpenMP for the GPW pair loops (PG_Cart_MnD evaluator).  Scoped to THIS lib's own TUs: -fopenmp on compile,
# the runtime on link, and QCHEM_OPENMP defined so the (opt-in, serial-by-default) parallel code activates.
# Prefer a detected OpenMP (libomp); else fall back to driving libgomp directly (this LLVM install lacks
# libomp).  This shared lib carries the runtime as a NEEDED entry, so consumers link cleanly.
if(OpenMP_CXX_FOUND)
    target_link_libraries(qcMolecule_BS PRIVATE OpenMP::OpenMP_CXX)
else()
    target_compile_options(qcMolecule_BS PRIVATE -fopenmp=libgomp)
    target_link_options   (qcMolecule_BS PRIVATE -fopenmp=libgomp)
endif()
target_compile_definitions(qcMolecule_BS PRIVATE QCHEM_OPENMP)
add_subdirectory(PG_Cart)
add_subdirectory(PG_Spherical)
add_subdirectory(PG_LibCint)
add_subdirectory(Evaluators)
add_subdirectory(tests)
add_subdirectory(bench)