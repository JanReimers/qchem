add_library(qcHamiltonian) 
target_sources(qcHamiltonian PRIVATE 
    Imp/Energy.C 
    Imp/Factory.C
    Internal/Imp/DiracKinetic.C
    Internal/Imp/RestMass.C 
    Internal/Imp/Ven.C
    Internal/Imp/Ven_Forces.C             # Ven's force (GeometryOptimizationPlan §13.2; 2026-11-08)
    Internal/Imp/Vee_FitDunlap_Forces.C   # the Dunlap Hartree force (GeometryOptimizationPlan §13.2; 2026-11-08)
    Internal/Imp/Vee_Fit_Forces.C         # the first-order fitted Hartree force (GeometryOptimizationPlan §13.2; 2026-11-09)
    Internal/Imp/PP_Forces.C              # the two PP terms' forces (GeometryOptimizationPlan §13.2; 2026-11-08)
    Internal/Imp/HF_Forces.C              # the HF Coulomb and exchange forces (GeometryOptimizationPlan §13.2; 2026-11-08)
    Internal/Imp/PP_Local.C
    Internal/Imp/PP_NonLocal.C
    Internal/Imp/HF_HT.C
    Internal/Imp/Vee.C
    Internal/Imp/Vxc.C
    Internal/Imp/Vee_Fit.C
    Internal/Imp/FittedVxc.C
    # The periodic (plane-wave / GPW) term set -- ONE module interface, five implementation units, the
    # same shape Internal/Terms.C has always had.  Split 2026-09-08 out of a single 1213-line TU.
    Internal/Imp/PWTerms_PP.C
    Internal/Imp/PWTerms_PP_Forces.C      # the PW PP structure-factor forces (GeometryOptimizationPlan §13.2; 2026-11-08)
    Internal/Imp/PWTerms_XC.C
    Internal/Imp/Vxc_Forces.C             # the XC term's force at fixed D (GeometryOptimizationPlan §13.2; lane B)
    Internal/Imp/Hubbard.C                # DFT+U (programme step 5, 2026-09-19)
    Internal/Imp/ACBN0.C                  # DFT+U increment 3: the ACBN0 (U,J) estimator (2026-09-21)
    # ⚠ THE SAMPLING ENGINE -- its OWN module since 2026-09-08 (qchem.Hamiltonian.Internal.DensitySampler).
    # It is not "H_ij and E from a density" and does not belong in this library; the module boundary is the
    # step that makes the library move a CMake change instead of surgery.  The HOME is an open decision:
    # it cannot go to qcFitting (that would close a cycle -- see the module header).
    # doc/CleanupCandidates.md R1.0e.
        Internal/Imp/HamiltonianImp.C
    Internal/Imp/HamiltonianForces.C      # the force fold over the terms (GeometryOptimizationPlan §13; 2026-11-08)
    Internal/Imp/Hamiltonians.C 
    Internal/Imp/SlaterExchange.C         # e_x, and Enzyme's v_x / f_x (the plugin rides this file)
    Internal/Imp/VWN_Correlation.C        # e_c, and Enzyme's v_c / f_c (the plugin rides this file)
    Internal/Imp/Libxc_LDA.C
    )
target_sources(qcHamiltonian PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Energy.C  
    Hamiltonian.C 
    HubbardEstimator.C
    HubbardChannels.C                     # the +U projectors, public (LinearResponsePlan R0)
    HubbardUTarget.C                      # the writes a U estimator makes to the +U term (LinearResponsePlan §3d Q6)
    TransitionFock.C                      # δF, the response kernel's answer (LinearResponsePlan C2)
    Factory.C
    Types.C
    Internal/HamiltonianImp.C 
    Internal/RealSpaceResidual.C
    Internal/GeometryGradient.C           # the term face of the analytic force (GeometryOptimizationPlan §13.2; 2026-11-08)
    Internal/HamiltonianTerm.C 
    Internal/Hamiltonians.C 
    Internal/ExchangeFunctional.C  
    Internal/SlaterExchange.C
    Internal/VWN_Correlation.C
    Internal/Hubbard.C                    # DFT+U: Hubbard_U + LowdinProjector (born on MatrixForward/Adjoint)
    Internal/ACBN0.C                      # DFT+U increment 3: ACBN0 realises HubbardUEstimator over HubbardProjection
    Internal/Libxc_LDA.C
    Internal/IonIon.C
    Internal/Kinetic.C
    Internal/Terms.C
      Internal/PWTerms.C
)

   
set_property(TARGET qcHamiltonian PROPERTY CXX_STANDARD 20)
target_include_directories(qcHamiltonian PRIVATE 
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_LIBXC_BUILD}
    )
target_link_libraries(qcHamiltonian PUBLIC qcBasisSet qcFitting qcOrbitals qcSymmetry qcStructure qcMesh qcMath qcCommon qcLASolver qcPseudopotential xc)
target_compile_definitions(qcHamiltonian PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})



add_subdirectory(tests)

# The LDA functionals' potentials and kernels are Enzyme's derivatives of their energies (Enzyme REQUIRED, family 2,
# 2026-11-09); only these two files carry the plugin.
qchem_enzyme_sources(qcHamiltonian Internal/Imp/SlaterExchange.C Internal/Imp/VWN_Correlation.C)
