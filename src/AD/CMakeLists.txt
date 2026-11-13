# qcAD: automatic differentiation by Enzyme (doc/GeometryOptimizationPlan.md, the AD spike; lane B, 2026-11-08).
#   - qchem.AD           : what a client sees -- NoGradient (the typed wrapper lives in the Imp TUs: see below)
#   - qchem.AD.XCKernel  : the Becke-mesh delta XC energy at fixed D as a PURE kernel of (R, D), and its AD gradients
#
# WHY ONE LIBRARY HOLDS THE KERNELS.  Enzyme differentiates LLVM IR, so it sees only bodies in the TU that asks for the derivative
# (or in one LTO link): a call into another TU, another shared library, or through a vtable is refused at COMPILE time (the spike's
# probes).  So a differentiable kernel and its __enzyme_autodiff call sit in one Imp TU here, and only this target carries the plugin.
add_library(qcAD)
target_sources(qcAD PRIVATE Imp/XCKernel.C)
target_sources(qcAD PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    AD.C
    XCKernel.C
)
set_property(TARGET qcAD PROPERTY CXX_STANDARD 20)
target_compile_definitions(qcAD PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_include_directories(qcAD PRIVATE ${QCHEM_BLAZE_DIR})
target_link_libraries(qcAD PUBLIC qcMath qcCommon)
# The plugin rides every qcAD TU, per file (qchem_enzyme_sources, root CMakeLists.txt: the plugin and its headers at the tree's codegen
# flags, and the plugin .so as an OBJECT_DEPENDS so a rebuilt plugin recompiles them).  The spike's per-target -O2
# -mprefer-vector-width=128 dodged EnzymeAD/Enzyme#3747 and #3748, both fixed at the pin cdd0d1c2 (2026-11-13).  Repros and the record:
# doc/Records/EnzymeSpike.md.  No flag set makes Enzyme safe: every AD gradient is gated by a finite-difference test.
qchem_enzyme_sources(qcAD Imp/XCKernel.C AD.C XCKernel.C)

add_subdirectory(tests)
