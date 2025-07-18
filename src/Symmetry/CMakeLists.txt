add_library(qcSymmetry) 
target_sources(qcSymmetry PRIVATE Atom_ECImp.C  ElectronCountsImp.C Irrep_QNsImp.C  Molecule_ECImp.C  Okmj.C  Orbital_QNsImp.C  SymmetryImp.C  UnitImp.C YlImp.C  YlmImp.C OkmjImp.C)
target_sources(qcSymmetry PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Spin.C 
    Symmetry.C 
    Irrep_QNs.C 
    Orbital_QNs.C 
    ElectronCounts.C 
    Angular.C 
    Yl.C Ylm.C Okmj.C
    Unit.C
    ElectronConfiguration.C
    Atom_EC.C 
    Molecule_EC.C 
    ) 
target_include_directories(qcSymmetry PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)
target_link_libraries(qcSymmetry PRIVATE qcCommon omlm)
set_property(TARGET qcSymmetry PROPERTY CXX_STANDARD 20)
# This linker association is required in order for cmake/gcc to find .gcm binaru module  interfact (BMI) files!!!!!
# target_link_libraries(qcMesh PUBLIC qcMesh qcCommon) 
  