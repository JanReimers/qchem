add_library(qcHamiltonian) 
target_sources(qcHamiltonian PRIVATE 
    Imp/Energy.C 
    Imp/Factory.C
    Internal/Imp/Kinetic.C  
    Internal/Imp/DiracKinetic.C 
    Internal/Imp/RestMass.C 
    Internal/Imp/Ven.C 
    Internal/Imp/Vee.C 
    Internal/Imp/Vxc.C 
    Internal/Imp/VxcPol.C 
    Internal/Imp/Vnn.C 
    Internal/Imp/FittedVxcPol.C 
    Internal/Imp/FittedVee.C 
    Internal/Imp/FittedVxc.C 
    Internal/Imp/HamiltonianTerm.C  
    Internal/Imp/HamiltonianImp.C  
    Internal/Imp/Hamiltonians.C 
    Internal/Imp/ExchangeFunctional.C  
    Internal/Imp/SlaterExchange.C 
    Internal/Imp/LDAVxc.C
    )
target_sources(qcHamiltonian PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Energy.C  
    Hamiltonian.C 
    Factory.C
    Internal/HamiltonianImp.C 
    Internal/HamiltonianTerm.C 
    Internal/Hamiltonians.C 
    Internal/ExchangeFunctional.C  
    Internal/SlaterExchange.C 
    Internal/LDAVxc.C 
    Internal/Terms.C
)

target_include_directories(qcHamiltonian PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Hamiltonian  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )
    
set_property(TARGET qcHamiltonian PROPERTY CXX_STANDARD 20)
target_link_libraries(qcHamiltonian PUBLIC qcChargeDensity qcFitting qcBasisSet qcSymmetry qcCluster qcMesh qcCommon)
