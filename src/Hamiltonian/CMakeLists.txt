add_library(qcHamiltonian) 
target_sources(qcHamiltonian PRIVATE 
    Imp/Energy.C 
    Internal/Imp/HamiltonianImp.C 
    Imp/HamiltonianTerm.C   
    Hamiltonians.C 
    Kinetic.C  
    DiracKinetic.C 
    RestMass.C 
    Ven.C 
    Vee.C 
    Vxc.C 
    VxcPol.C 
    FittedVxcPol.C 
    FittedVee.C 
    FittedVxc.C 
    LDAVxc.C Vnn.C 
    Internal/Imp/ExchangeFunctional.C  
    SlaterExchange.C 
    ExchangeFunctionalFactories.C 
    Factory.C
    )
target_sources(qcHamiltonian PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Energy.C  
    Hamiltonian.C 
    Internal/HamiltonianImp.C 
    Internal/ExchangeFunctional.C  
)

target_include_directories(qcHamiltonian PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Hamiltonian  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )
    
set_property(TARGET qcHamiltonian PROPERTY CXX_STANDARD 20)
target_link_libraries(qcHamiltonian PUBLIC qcChargeDensity qcFitting qcBasisSet qcSymmetry qcCluster qcMesh qcCommon)
