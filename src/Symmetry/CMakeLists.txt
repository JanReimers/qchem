add_library(qcSymmetry) 
target_sources(qcSymmetry PRIVATE 
    Imp/Atom_EC.C  
    Imp/ElectronCounts.C 
    Imp/Irrep_QNs.C  
    Imp/Molecule_EC.C  
    Imp/Orbital_QNs.C  
    Imp/Symmetry.C  
    Imp/Unit.C 
    Imp/Yl.C  
    Imp/Ylm.C 
    Imp/Okmj.C
    Imp/BlochQN.C
    )
target_sources(qcSymmetry PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Spin.C 
    Symmetry.C 
    Irrep_QNs.C 
    Orbital_QNs.C 
    ElectronCounts.C 
    Angular.C 
    Yl.C Ylm.C Okmj.C BlochQN.C
    Unit.C
    ElectronConfiguration.C
    Atom_EC.C 
    Molecule_EC.C 
    ) 
# This linker association is required in order for cmake/gcc to find .gcm binaru module  interfact (BMI) files!!!!!
target_link_libraries(qcSymmetry PRIVATE qcCommon omlm)
set_property(TARGET qcSymmetry PROPERTY CXX_STANDARD 20)
  