// File: Response/OrbitalFrame.C  The AO <-> MO bridge of a linear response (doc/LinearResponsePlan.md §3c).
//
// WHY A SEPARATE OBJECT FROM THE REFERENCE.  The Reference is R0 in the ORBITAL basis -- scalar-agnostic
// (complex throughout), validated in R0, and it never needed a coefficient.  The response KERNEL speaks the
// AO basis of the Hamiltonian's terms, which is T-typed (real for a molecule, complex for a general k-point).
// This frame is the one place that holds each block's coefficients C and basis, and it answers the two
// OPERATIONS that cross the boundary (the IrrepCD preference: no GetC()):
//   ToAO: δD in the MO basis  ->  a TransitionDensity the kernel can consume,  δD_AO = C δD_MO C^†
//   ToMO: a TransitionFock    ->  δF in the MO basis R0 consumes,             δF_MO = C^† δF_AO C
// It is built block-aligned with its Reference (same order, same orbital counts -- checked).
//
// q = 0 ONLY in R1: every block is its own partner, so δD_AO is square, and Hermitian for a Hermitian
// perturbation.  On a REAL block (T = double) the MO δD must be real: its imaginary part is checked, not
// dropped silently (the real-TRIM rule, doc/RealComplexPlan.md).
module;
#include <memory>
#include <vector>
export module qchem.Response.OrbitalFrame;
export import qchem.Response.Reference;
export import qchem.ChargeDensity.TransitionDensity;
export import qchem.Hamiltonian.TransitionFock;

export namespace qchem::Response
{

//! One block of the frame: the block's irrep, its AO basis and its orbital coefficients (AO x orbital).
template <class T> struct FrameBlock
{
    Irrep                            irrep;
    const ChargeDensity::tobs_t<T>*  bs=nullptr;
    mat_t<T>                         C;
};

template <class T> class OrbitalFrame
{
public:
    //! THROWS unless \a blocks line up with \a ref's (irrep order and orbital count per block).
    OrbitalFrame(const Reference& ref, std::vector<FrameBlock<T>> blocks);
    //! δD (MO basis, on \a rule's pairs) -> the AO transition density.  THROWS if a block is not its own
    //! partner (q != 0 is R3's), or a real block receives a complex δD.
    std::unique_ptr<ChargeDensity::TransitionDensity<T>> ToAO(const BlockPairs& dD,
                                                              std::shared_ptr<const Symmetry::SelectionRule> rule) const;
    //! δF (AO basis) -> MO block pairs on \a rule's pairs: \f$C_{\rm bra}^\dagger\,\delta F\,C_{\rm ket}\f$.
    BlockPairs ToMO(const Hamiltonian::TransitionFock<T>& dF, const Symmetry::SelectionRule& rule) const;
    size_t NumBlocks() const {return itsBlocks.size();}
private:
    const Reference&             itsRef;
    std::vector<FrameBlock<T>>   itsBlocks;
};

} // namespace
