// File: AD/AD.C  Automatic differentiation by Enzyme: what a client of a geometry gradient sees (lane B, 2026-11-08).
//
// Enzyme (an LLVM plugin) differentiates a function's IR at compile time.  It is REQUIRED (user ruling 2026-11-09; the root
// CMakeLists.txt stops configure without it), so an AD derivative never refuses for want of it.  No hand-coded derivative is kept
// beside an AD one as a fallback.
//
// THE TYPED WRAPPER is not exported: Enzyme's entry point takes the function as an untyped pointer (its own header hides that behind
// enzyme::autodiff), and the plugin must compile the TU that instantiates it, so the wrapper lives in qcAD's Imp TUs beside the kernels
// it differentiates (src/AD/Imp/XCKernel.C).  A client asks a KERNEL module (qchem.AD.XCKernel) for its gradient.
module;
#include <string>
export module qchem.AD;

export namespace qchem::AD
{

//! Why a geometry gradient is not available: the refusal of a derivative face whose basis or term has none (an Outcome's error).
//! \c why reads as a sentence ("LAPW gives no overlap pullback yet").
struct NoGradient { std::string why; };

} // namespace qchem::AD
