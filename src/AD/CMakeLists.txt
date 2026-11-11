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
# PRIVATE: the plugin changes only how qcAD's own TUs compile; its interface (and so every importer's .pcm) does not see it.
# -O2 -mprefer-vector-width=128 (after the tree's -O3, so they win): Enzyme 1a87ccc on LLVM 21.1.6 SEGFAULTED at -O2 -march=native
# (256-bit vectors, AVX2; EnzymeAD/Enzyme#3748, open at the pin) on the spike's blocked kernel, and its loop-cache defect gave a silently
# wrong dE/dD at -O3 there (fixed by v0.0.308; #3747 is a second such shape, open).  -march itself must stay (a .pcm and its importers
# agree on target features).  Repros and the record: doc/Records/EnzymeSpike.md.  No flag set makes Enzyme safe: every AD gradient is
# gated by a finite-difference test.
target_compile_options(qcAD PRIVATE "-fplugin=${QCHEM_ENZYME_PLUGIN}" -O2 -mprefer-vector-width=128)
target_include_directories(qcAD PRIVATE "${QCHEM_ENZYME_INCLUDE}")

add_subdirectory(tests)
