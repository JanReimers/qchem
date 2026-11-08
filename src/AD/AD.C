// File: AD/AD.C  Automatic differentiation by Enzyme: what a client of an AD gradient sees (lane B, 2026-11-08).
//
// Enzyme (an LLVM plugin) differentiates a function's IR at compile time.  It is TOOLCHAIN and OPTIONAL (controller ruling 2026-11-08):
// a build without it (QCHEM_ENZYME=OFF) still builds every kernel, and every AD gradient REFUSES with NoGradient, naming why.  No
// hand-coded derivative is kept beside an AD one as a fallback.
//
// THE TYPED WRAPPER is not exported: Enzyme's entry point takes the function as an untyped pointer (its own header hides that behind
// enzyme::autodiff), and the plugin must compile the TU that instantiates it, so the wrapper lives in qcAD's Imp TUs beside the kernels
// it differentiates (src/AD/Imp/AD.C).  A client asks a KERNEL module (qchem.AD.XCKernel) for its gradient, as an Outcome.
module;
#include <string>
export module qchem.AD;

export namespace qchem::AD
{

//! Why an AD gradient is not available.  \c why reads as a sentence ("this build has no Enzyme (QCHEM_ENZYME=OFF)").
struct NoGradient { std::string why; };

//! Was this build compiled with Enzyme (QCHEM_ENZYME=ON)?  For a TEST that checks the refusal or the numbers; a client never branches
//! on it -- it reads the Outcome.
bool Available();

} // namespace qchem::AD
