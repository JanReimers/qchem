add_library(qcSymmetry)
target_sources(qcSymmetry PRIVATE
    # --- root: shared sym_t spine ---
    Imp/Symmetry.C
    Imp/Irrep.C
    Imp/Orbital_QNs.C
    Imp/Factory.C
    Imp/Unit.C
    # --- Atom/ ---
    Atom/Imp/Spherical.C
    Atom/Internal/Imp/SphericalQNs.C   # merged: was Yl.C + Ylm.C + Okmj.C (all impl units of one module)
    # --- Molecule/ ---
    Molecule/Imp/Irrep.C
    Molecule/Imp/PointGroup.C
    Molecule/Imp/CartesianRep.C
    Molecule/Imp/SphericalRep.C
    Molecule/Imp/OperationRep.C
    Molecule/Imp/CharacterTable.C
    Molecule/Imp/AbelianGroup.C
    Molecule/Imp/SALC.C
    # --- Lattice_3D/ ---
    Lattice_3D/Imp/BlochQN.C
    Lattice_3D/Imp/SpaceGroup.C
    Lattice_3D/Imp/BZReduction.C
    )
target_sources(qcSymmetry PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    # --- root: shared sym_t spine ---
    Spin.C
    Symmetry.C
    Irrep.C
    Orbital_QNs.C
    Factory.C
    Unit.C
    # --- Atom/ ---
    Atom/Spherical.C
    Atom/Internal/SphericalQNs.C
    # --- Molecule/ ---
    Molecule/Irrep.C
    Molecule/PointGroup.C
    Molecule/ShellRep.C
    Molecule/CartesianRep.C
    Molecule/SphericalRep.C
    Molecule/OperationRep.C
    Molecule/CharacterTable.C
    Molecule/AbelianGroup.C
    Molecule/SALC.C
    # --- Lattice_3D/ ---
    Lattice_3D/BlochQN.C
    Lattice_3D/SpaceGroup.C
    Lattice_3D/BZReduction.C
    )

target_compile_options(qcSymmetry PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcSymmetry PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_include_directories(qcSymmetry PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/submodules/blaze)
# This linker association is required in order for cmake/gcc to find .gcm binaru module  interfact (BMI) files!!!!!
target_link_libraries(qcSymmetry PUBLIC qcMath qcCommon)
set_property(TARGET qcSymmetry PROPERTY CXX_STANDARD 20)

add_subdirectory(tests)