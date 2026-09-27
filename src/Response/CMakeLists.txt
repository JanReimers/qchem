# qcResponse -- linear response (DFPT / CPHF) as GENERAL perturbation theory (doc/LinearResponsePlan.md).
# Sits ABOVE qcWaveFunction and qcHamiltonian (it reads a converged state and the Hamiltonian's capabilities)
# and BELOW the qchem aggregate (the SolidCalculation facade drives it).  Stage R0: the independent-particle
# response chi0 by sum over states, over the +U projectors.
add_library(qcResponse)
target_sources(qcResponse PRIVATE Imp/Reference.C Imp/Probe.C Imp/Response.C)
target_sources(qcResponse PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Reference.C
    Probe.C
    Response.C
)
target_include_directories(qcResponse PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    )
set_property(TARGET qcResponse PROPERTY CXX_STANDARD 20)
target_link_libraries(qcResponse PUBLIC qcWaveFunction qcOrbitals qcHamiltonian qcBasisSet qcElConfig qcSymmetry qcMath qcCommon)
target_compile_definitions(qcResponse PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})

add_subdirectory(tests)
