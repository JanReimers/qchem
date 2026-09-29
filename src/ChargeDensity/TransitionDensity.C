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

//! One block PAIR of an AO transition density: the (k+q, σ) bra block, the (k, σ) ket block, and δD on the pair
//! (bra rows x ket columns).  At q = 0 (or a totally symmetric perturbation) bra == ket, and δD is square --
//! Hermitian for a Hermitian perturbation.  At q != 0 it is neither square in meaning nor Hermitian: its conjugate
//! partner lives on the (k, k+q) pair (doc/LinearResponsePlan.md §3d, C1).  A construction VALUE.
template <class T> struct TransitionBlock
{
    Irrep            bra, ket;
    const tobs_t<T>* braBs=nullptr;
    const tobs_t<T>* ketBs=nullptr;
    mat_t<T>         dD;
};

//! \brief THE HARTREE FACE of a PERIODIC transition density (doc/LinearResponsePlan.md §3d, C1): \f$\delta V_H(G+q)\f$
//! summed over its block pairs, on the CD fit basis's ball.  A cross-cast capability, like \c FourierDensity is on a
//! ground-state density -- and a DIFFERENT type from it on purpose: its answer is a \c ΔGq_Map (ruling Q7), and it
//! is never star-averaged (a perturbation breaks the imposed group, §3d finding 5).
class TransitionFourierDensity
{
public:
    virtual ~TransitionFourierDensity() = default;
    virtual ΔGq_Map GetTransitionRepulsion(const BasisSet::cFIT_CD_ABS&) const=0;
};

//! \brief What an XC sampler projects a periodic transition density THROUGH -- the pair sibling of
//! \c Fitting::ScalarProjector.  The sampler implements it over its own quadrature (Φ tables on the δ/Becke mesh, or
//! the raw collocation on a raster); the transition density contracts each of its own pairs into it
//! (\c ProjectableTransition), so δD never leaves the density.  \a q is the density's one wave vector.
class TransitionProjector
{
public:
    virtual ~TransitionProjector() = default;
    //! This quadrature's values of ONE pair's transition density (bra x ket δD), one per coefficient.
    virtual cvec_t Forward(const tobs_t<dcmplx>& bra, const tobs_t<dcmplx>& ket, const mat_t<dcmplx>& dD,
                           const rvec3_t& q) const=0;
    virtual size_t NumCoefficients() const=0;   //!< the length of every vector \c Forward returns
};

//! \brief "I can be projected pair by pair" -- the periodic transition density's face onto a \c TransitionProjector:
//! \f$\sum_{\rm pairs}\f$ \c Forward(bra, ket, δD, q).  NEVER symmetrized (§3d finding 5).
class ProjectableTransition
{
public:
    virtual ~ProjectableTransition() = default;
    virtual cvec_t ProjectOnto(const TransitionProjector&) const=0;
};

//! The wave vector a selection rule shifts by (fractional; its \c Symmetry::WaveVectorShift face), or 0 for a rule
//! without one -- the ONE place that question is answered, for the transition density and the terms alike.
inline rvec3_t WaveVectorOf(const Symmetry::SelectionRule& rule)
{
    auto* w=dynamic_cast<const Symmetry::WaveVectorShift*>(&rule);   // abstract -> abstract
    return w ? w->q() : rvec3_t(0,0,0);
}

//! \brief The AO-matrix transition density over \a blocks, coupled by \a rule.  THROWS unless \a rule couples every
//! block's ket to its bra.  The wave vector comes from the rule (\c Symmetry::WaveVectorShift; none = q = 0).
//! On the real path (T = double, q = 0 only) the result carries the HF sweep face \c tHF_System_CD<double>; on the
//! periodic path (T = dcmplx) it carries \c TransitionFourierDensity and \c ProjectableTransition.
template <class T> std::unique_ptr<TransitionDensity<T>>
AO_TransitionDensity_Factory(std::vector<TransitionBlock<T>> blocks, std::shared_ptr<const Symmetry::SelectionRule> rule);

using rTransitionDensity=TransitionDensity<double>;   using cTransitionDensity=TransitionDensity<dcmplx>;

} // namespace
