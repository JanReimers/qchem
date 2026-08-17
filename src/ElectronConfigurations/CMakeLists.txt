add_library(qcElConfig) 
target_sources(qcElConfig PRIVATE
    Internal/Imp/ElectronCounts.C
    Imp/Atom_EC.C
    Imp/AtomDirac_EC.C
    Imp/Molecule_EC.C
    Imp/MolecularSym_EC.C
    Imp/Crystal_EC.C
    Imp/OccupationPolicy.C
    )
target_sources(qcElConfig PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Internal/ElectronCounts.C
    ElectronConfiguration.C
    OrbitalView.C
    OccupationPolicy.C
    Atom_EC.C
    AtomDirac_EC.C
    Molecule_EC.C
    MolecularSym_EC.C
    Crystal_EC.C
    )

target_compile_options(qcElConfig PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcElConfig PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_include_directories(qcElConfig PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/submodules/blaze)
# This linker association is required in order for cmake/gcc to find .gcm binaru module  interfact (BMI) files!!!!!
# The OccupationPolicy reads orbitals through its own OrbitalView DIP face (qcOrbitals sits ABOVE this
# library and implements it -- V1.11 inc 3, user ruling: invert, never re-point the DAG).  qcMath supplies
# the coefficient/reference types (vec_t/mat_t via qchem.Types).
target_link_libraries(qcElConfig PUBLIC qcSymmetry qcMath)
set_property(TARGET qcElConfig PROPERTY CXX_STANDARD 20)

add_subdirectory(tests)