add_library(qcSymmetry) 
target_sources(qcSymmetry PRIVATE 
    Imp/Spherical.C
    Imp/Irrep.C  
    Imp/Orbital_QNs.C  
    Imp/Symmetry.C
    Imp/Factory.C
    Imp/PointGroup.C
    Imp/CartesianRep.C
    Imp/CharacterTable.C
    Imp/AbelianGroup.C
    Imp/SALC.C
    Internal/Imp/Unit.C
    Internal/Imp/Yl.C  
    Internal/Imp/Ylm.C 
    Internal/Imp/Okmj.C
    Internal/Imp/BlochQN.C
    )
target_sources(qcSymmetry PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Spin.C 
    Symmetry.C 
    Irrep.C 
    Orbital_QNs.C 
    Spherical.C 
    Factory.C
    PointGroup.C
    CartesianRep.C
    CharacterTable.C
    AbelianGroup.C
    SALC.C
    Internal/Spherical.C
    Internal/BlochQN.C
    Internal/Unit.C

    ) 

target_compile_options(qcSymmetry PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcSymmetry PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_include_directories(qcSymmetry PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/submodules/blaze)
# This linker association is required in order for cmake/gcc to find .gcm binaru module  interfact (BMI) files!!!!!
target_link_libraries(qcSymmetry PUBLIC qcCommon)
set_property(TARGET qcSymmetry PROPERTY CXX_STANDARD 20)
  