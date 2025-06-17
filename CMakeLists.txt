cmake_minimum_required(VERSION 3.5.0)
project(qchem VERSION 0.1.0 LANGUAGES CXX)

enable_language (Fortran)

add_library(qchem) 
target_sources(qchem PRIVATE src/SCFIterator.C)
add_subdirectory(./src)

# -Wno-unknown-pragmas
target_compile_options(qchem PRIVATE -Wall -Wno-unknown-pragmas
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_compile_definitions(qchem PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_include_directories(qchem PRIVATE include ./submodules/oml/include )
target_link_libraries(qchem qcWaveFunction qcSCFAccelerator qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcMolecule_BS qcAtom_BS qcBasisSet qcLASolver qcCluster qcMesh qcSymmetry  qcCommon  oml lapack gtest wignerSymbols)
add_subdirectory(./UnitTests)
add_subdirectory(./gtkapp)
add_subdirectory(./submodules/oml)
add_subdirectory(./submodules/wignerSymbols)



