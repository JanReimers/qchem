add_library(qcSymmetry) 
target_sources(qcSymmetry PRIVATE 
    Imp/Spherical.C
    Imp/Atom_EC.C  
    Imp/ElectronCounts.C 
    Imp/Irrep_QNs.C  
    Imp/Molecule_EC.C  
    Imp/Orbital_QNs.C  
    Imp/Symmetry.C  
    Imp/Factory.C
    Internal/Imp/Unit.C 
    Internal/Imp/Yl.C  
    Internal/Imp/Ylm.C 
    Internal/Imp/Okmj.C
    Internal/Imp/BlochQN.C
    )
target_sources(qcSymmetry PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Spin.C 
    Symmetry.C 
    Irrep_QNs.C 
    Orbital_QNs.C 
    ElectronCounts.C 
    Spherical.C 
    ElectronConfiguration.C
    Atom_EC.C 
    Atom_Dirac_EC.C
    Molecule_EC.C 
    Factory.C
    Internal/Spherical.C 
    Internal/BlochQN.C
    Internal/Unit.C

    ) 

target_compile_options(qcSymmetry PRIVATE -Wall 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -pg>
)
target_include_directories(qcSymmetry PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/submodules/blaze)
# This linker association is required in order for cmake/gcc to find .gcm binaru module  interfact (BMI) files!!!!!
target_link_libraries(qcSymmetry PRIVATE qcCommon)
set_property(TARGET qcSymmetry PROPERTY CXX_STANDARD 20)
  