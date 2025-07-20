add_library(qcOrbitals) 
target_sources(qcOrbitals PRIVATE Internal/Imp/Orbital.C  Internal/Imp/TOrbital.C  Internal/Imp/TOrbitals.C Imp/Factory.C)
target_sources(qcOrbitals PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Orbitals.C 
    Factory.C
    Internal/TOrbital.C  
    Internal/TOrbitals.C
)
target_include_directories(qcOrbitals PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Orbitals  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )
set_property(TARGET qcOrbitals PROPERTY CXX_STANDARD 20)
target_link_libraries(qcOrbitals PUBLIC qcChargeDensity qcBasisSet qcSymmetry qcCluster qcMesh qcCommon)
