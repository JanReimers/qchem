// File: Symmetry/Imposition.C  What a run IMPOSES, as its consumers ask it, and the spin-pair projector (V-BS-SYMMETRY S3).
//
// The interface is structure-neutral on purpose (user ruling B, 2026-11-03): a crystal's imposition realizes it today
// (qchem::CrystalImposition, resolved from the lattice and the options), and a molecule's imposed point group can realize it
// later (V-FACADE-PARITY) -- a run without one hands nullptr, which is a stage, not a design.
module;
#include <cassert>
#include <memory>
#include <vector>
export module qchem.Symmetry.Imposition;
export import qchem.Symmetry.Representation;
export import qchem.Symmetry.StarGroup;   // StarGroup (the density star average)
import qchem.Blaze;

export namespace qchem::Symmetry
{

//! \brief The run's imposed symmetry, answered for each consumer.  Data-free interface.
class SymmetryImposition
{
public:
    virtual ~SymmetryImposition() = default;
    //! The σ-carrying group the (D↑, D↓) PAIR is projected under ({} = no magnetic pair projection).  Only ops that map every
    //! orbital block onto ITSELF belong here: a flip op relates the two channels, so only a caller holding both can impose it.
    virtual const std::vector<SymOp>& SpinPairOps() const=0;
    //! The group the composite density star-averages under (ρ̃ in G space and the raster ρ alike), or null = the trivial {E}.
    //! Shared, because every density built under this imposition holds it.
    virtual std::shared_ptr<const StarGroup> DensityGroup() const=0;
};

//! \brief Impose \a imp's spin-pair group on a block's (D↑, D↓), in place:
//! \f$D_s \leftarrow |G|^{-1}\sum_g \mathrm{AdjointImage}(g, D_{\sigma_g(s)})\f$, a flip op reading the OTHER channel; then
//! the Hermitian part (one rounding, at the end).  A projector because the set is a closed group (the imposition gate) and
//! g -> M_g a representation on the block.  It commutes with every Hamiltonian term built from D, so imposing it here
//! imposes it on ρ, m, the mixer's history and the +U occupations at once.  A no-op when \a imp has no spin-pair group.
template <class T> void ProjectSpinPair(const SymmetryImposition& imp, const Representation<T>& rep, hmat_t<T>& Dup, hmat_t<T>& Ddn)
{
    const std::vector<SymOp>& ops=imp.SpinPairOps();
    if (ops.empty()) return;
    assert(Ddn.rows()==Dup.rows());
    const size_t n=Dup.rows();
    mat_t<T> up(n,n,T(0)), dn(n,n,T(0));
    for (const SymOp& g : ops)
    {
        const bool flip=g.sigma==SpinAction::Flip;
        rep.AddAdjointImage(g, flip ? Ddn : Dup, up);
        rep.AddAdjointImage(g, flip ? Dup : Ddn, dn);
    }
    const double w=1.0/double(ops.size());
    // Hermitian to roundoff by construction; store the Hermitian part (an hmat_t refuses a 1e-16 asymmetry).
    Dup=blazem::declherm(0.5*w*(up+blazem::ctrans(up)));
    Ddn=blazem::declherm(0.5*w*(dn+blazem::ctrans(dn)));
}

} // namespace qchem::Symmetry
