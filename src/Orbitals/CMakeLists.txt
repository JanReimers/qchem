add_library(qcOrbitals) 
target_sources(qcOrbitals PRIVATE Orbital.C  TOrbital.C  TOrbitals.C)
target_include_directories(qcOrbitals PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Orbitals  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )