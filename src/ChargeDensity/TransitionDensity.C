// File: ChargeDensity/TransitionDensity.C  The first-order density of a linear response
// (doc/LinearResponsePlan.md §3 row C1, §3c).
//
// WHAT IT IS.  δD, the change of the density matrix under a perturbation of symmetry Γ_pert: on the block
// pairs the perturbation couples (k+q <- k for a lattice wave vector; a block with itself at q = 0 or for a
// totally symmetric perturbation).  The response kernel (tHamiltonian::MakeResponseKernel) consumes it and
// answers the first-order Fock change δF = (∂F/∂D) δD.
//
// ★ IT IS NOT A tChargeDensity, AND NOT A tDM_CD (LSP; user question 2026-09-28).  δD is traceless and
// indefinite, is not normalisable and is not mixable along an SCF trajectory.  As a tDM_CD,
// tHamiltonian::GetMatrix(bs,s,δD) would COMPILE -- and for XC be silently wrong (the functionals skip
// ρ <= 0, dropping the negative half of δρ).  At q != 0 δD lives on a block PAIR and is no single irrep
// block at all, and a plane-wave Sternheimer code has no δD matrix, so the inheritance would force a second
// type at R3 anyway.  What a transition density IS, is narrower: at q = 0 it is a matrix that can be
// SCATTERED THROUGH THE ERI (tHF_System_CD) -- exactly what a LINEAR operator (J, K) consumes -- and it
// inherits that face and nothing more.
//
// ITS CAPABILITIES ARE CROSS-CAST FACES, EACH ARRIVING WITH THE STAGE THAT CONSUMES IT: the HF sweep
// (tHF_System_CD<double>, R1), δρ_σ on the XC mesh (R2), δρ(G+q) for Hartree (R3).  The REPRESENTATION is
// the concrete's business: AO_TransitionDensity holds AO block matrices (Gaussian, dense PW); a large-cutoff
// PW code would hold {ψ_v, δψ_v} instead (§4b).
module;
#include <memory>
#include <vector>
export module qchem.ChargeDensity.TransitionDensity;
export import qchem.Symmetry.SelectionRule;
export import qchem.Symmetry.Irrep;
export import qchem.Symmetry.Spin;
export import qchem.ChargeDensity.Types;   // tobs_t
import qchem.Types;

export namespace qchem::ChargeDensity
{

template <class T> class TransitionDensity
{
public:
    virtual ~TransitionDensity() = default;
    //! Which bra block each ket block couples to (the perturbation's symmetry label, S1).  Shared, because the
    //! Fock change it induces carries the same label.
    virtual std::shared_ptr<const Symmetry::SelectionRule> Coupling() const = 0;
    //! A transient freshness serial from the SAME clock as every density (NextDensityVersion): a term keys
    //! its response cache on it.  Distinct transition densities have distinct serials.
    virtual size_t Version() const = 0;
    //! The \a s channel as a VIEW (non-owning, alive as long as this); \c Spin::None answers this whole
    //! density.  NULL when \a s is not resolved -- the \c tSpinResolved_CD probe idiom.
    virtual const TransitionDensity* Channel(const Spin& s) const = 0;
};

//! One block of an AO transition density: the ket block's irrep and basis, and δD on it.  q = 0 only (a
//! block couples to ITSELF, so δD is square and, for a Hermitian perturbation, Hermitian); the (k+q, k) pair
//! form is R3's.  A construction VALUE.
template <class T> struct TransitionBlock
{
    Irrep            irrep;
    const tobs_t<T>* bs=nullptr;
    hmat_t<T>        dD;
};

//! \brief The AO-matrix transition density over \a blocks, coupled by \a rule.  THROWS unless \a rule couples
//! every block to itself (the q = 0 / totally-symmetric case; a q != 0 transition density is stage R3's).
//! On the real path (T = double) the result carries the HF sweep face \c tHF_System_CD<double>.
template <class T> std::unique_ptr<TransitionDensity<T>>
AO_TransitionDensity_Factory(std::vector<TransitionBlock<T>> blocks, std::shared_ptr<const Symmetry::SelectionRule> rule);

using rTransitionDensity=TransitionDensity<double>;   using cTransitionDensity=TransitionDensity<dcmplx>;

} // namespace
