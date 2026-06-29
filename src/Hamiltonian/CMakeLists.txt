add_library(qcHamiltonian) 
target_sources(qcHamiltonian PRIVATE 
    Imp/Energy.C 
    Imp/Factory.C
    Internal/Imp/Kinetic.C  
    Internal/Imp/DiracKinetic.C 
    Internal/Imp/RestMass.C 
    Internal/Imp/Ven.C
    Internal/Imp/PP_Local.C
    Internal/Imp/PP_NonLocal.C
    Internal/Imp/Vee.C
    Internal/Imp/Vxc.C 
    Internal/Imp/VxcPol.C 
    Internal/Imp/Vnn.C 
    Internal/Imp/FittedVxcPol.C 
    Internal/Imp/FittedVee.C 
    Internal/Imp/FittedVxc.C
    Internal/Imp/PWTerms.C
    Internal/Imp/HamiltonianImp.C
    Internal/Imp/Hamiltonians.C 
    Internal/Imp/ExchangeFunctional.C  
    Internal/Imp/SlaterExchange.C 
    Internal/Imp/Libxc_LDA_Exchange.C 
    Internal/Imp/LDAVxc.C
    )
target_sources(qcHamiltonian PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Energy.C  
    Hamiltonian.C 
    Factory.C
    Types.C
    Internal/HamiltonianImp.C 
    Internal/HamiltonianTerm.C 
    Internal/Hamiltonians.C 
    Internal/ExchangeFunctional.C  
    Internal/SlaterExchange.C
    Internal/VWN_Correlation.C
    Internal/Libxc_LDA_Exchange.C
    Internal/LDAVxc.C
    Internal/Terms.C
    Internal/PWTerms.C
)

   
set_property(TARGET qcHamiltonian PROPERTY CXX_STANDARD 20)
target_include_directories(qcHamiltonian PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_BINARY_DIR}/submodules/libxc
    )
target_link_libraries(qcHamiltonian PUBLIC qcBasisSet qcFitting qcOrbitals qcSymmetry qcStructure qcMesh qcMath qcCommon qcLASolver qcPseudopotential xc)
target_compile_options(qcHamiltonian PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcHamiltonian PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})


