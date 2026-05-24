add_library(qcOrbitals) 
target_sources(qcOrbitals PRIVATE Internal/Imp/Orbital.C  Internal/Imp/TOrbital.C  Internal/Imp/TOrbitals.C Imp/Factory.C Imp/EnergyLevel.C)
target_sources(qcOrbitals PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Orbitals.C 
    Factory.C
    EnergyLevel.C
    Types.C
    Internal/TOrbital.C  
    Internal/TOrbitals.C
)
set_property(TARGET qcOrbitals PROPERTY CXX_STANDARD 20)
target_include_directories(qcOrbitals PRIVATE 
      ${PROJECT_SOURCE_DIR}/submodules/blaze
)
target_link_libraries(qcOrbitals PUBLIC qcChargeDensity qcBasisSet1 qcSymmetry qcCluster qcMesh qcCommon)
target_compile_options(qcOrbitals PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -pg>
)
target_compile_definitions(qcOrbitals PUBLIC ${CXX_COMMON_DEFINES}
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
