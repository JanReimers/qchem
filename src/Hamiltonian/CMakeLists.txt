add_library(qcHamiltonian) 
target_sources(qcHamiltonian PRIVATE FittedVxcPol.C LDAVxc.C Vnn.C DiracKinetic.C Hamiltonian.C RestMass.C Vxc.C 
ExchangeFunctionalImp.C  Hamiltonians.C SlaterExchange.C VxcPol.C FittedVee.C HamiltonianTerm.C  Vee.C 
FittedVxc.C Kinetic.C  Ven.C ExchangeFunctionalFactories.C TotalEnergy.C Factory.C)
target_include_directories(qcHamiltonian PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/Hamiltonian  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )