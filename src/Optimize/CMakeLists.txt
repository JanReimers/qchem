# qcOptimize: minimise any cost with any algorithm (V-OPTIMIZE-DIP, doc/Records/VOptimizeDIP_Design.md; user ruling 2026-11-09).
#   - qchem.Optimize.Cost     : the Cost face and its capability faces (Projected, ModelHessian, ValueNoise), the typed tCost<Fail>
#   - qchem.Optimize.Minimize : the ask/tell Stepper, the client's StoppingRule, the Minimize driver
#   - qchem.Optimize.BFGS_RFO : BFGS + RFO + trust radius (moved from qchem.GeometryOptimizer)
#   - qchem.Optimize          : the umbrella
# NO CHEMISTRY: it links qcMath (Blaze, rvec_t) and qcCommon (Outcome, Reporting) only (user ruling 2026-11-09, decision 2).  Every
# library with something to minimise depends on it, never the converse; the linker enforces that.
add_library(qcOptimize)
target_sources(qcOptimize PRIVATE Imp/Minimize.C Imp/BFGS_RFO.C)
target_sources(qcOptimize PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Cost.C
    Minimize.C
    BFGS_RFO.C
    Optimize.C
)
set_property(TARGET qcOptimize PROPERTY CXX_STANDARD 20)
target_include_directories(qcOptimize PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_compile_definitions(qcOptimize PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qcOptimize PUBLIC qcMath qcCommon)

add_subdirectory(tests)
