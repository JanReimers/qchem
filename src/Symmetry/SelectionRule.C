// File: Symmetry/SelectionRule.C  Which symmetry block a perturbation couples to which
// (doc/LinearResponsePlan.md §3 row S1, §3c).
//
// A first-order perturbation of symmetry Γ_pert couples a ket block of symmetry Γ_ket only to the bra blocks
// in Γ_pert ⊗ Γ_ket.  That ONE rule covers every pairing linear response needs, so no perturbation invents
// its own:
//   * a lattice wave vector q couples k to k+q (q is the perturbation's irrep of the translation group) --
//     the concrete is Lattice_3D::MeshShift;
//   * a totally symmetric perturbation, or q = 0, couples every block to ITSELF -- Invariant, below;
//   * a molecular dipole in C2v couples A1 to B1 (a point-group product table: a later concrete);
//   * a spin flip is ΔM_s = ±1 (the SOC increment).
// Spin is NOT decided here (yet): the response Reference pairs same-spin blocks.
module;
export module qchem.Symmetry.SelectionRule;
export import qchem.Symmetry;

export namespace qchem::Symmetry
{

//! \brief Does a perturbation couple the \a ket block to the \a bra block?  The perturbation's symmetry label,
//! seen from the one question its clients ask of it.
class SelectionRule
{
public:
    virtual ~SelectionRule() = default;
    virtual bool Couples(const Symmetry& bra, const Symmetry& ket) const = 0;
};

//! \brief A perturbation that transforms as the totally symmetric irrep (or a q = 0 lattice perturbation):
//! every block couples to itself and to nothing else.  "Same block" is the same \c SequenceIndex, which is
//! the block identity within one run's symmetry.
class Invariant : public virtual SelectionRule
{
public:
    virtual bool Couples(const Symmetry& bra, const Symmetry& ket) const override
    {return bra.SequenceIndex()==ket.SequenceIndex();}
};

//! \brief A selection rule that IS a lattice wave-vector shift \f$k\to k+q\f$ -- a cross-cast CAPABILITY on the
//! rule (the RealBlock idiom), asked by the one client that needs the number: a periodic transition density,
//! which Fourier-transforms at \f$G+q\f$ and collocates its periodic part \f$e^{-iq\cdot r}\delta\rho\f$
//! (doc/LinearResponsePlan.md §3d B1/C1).  A rule WITHOUT it (\c Invariant, a point-group product) is a q = 0
//! perturbation as far as a lattice is concerned.  Its own face so \c SelectionRule stays theory- and
//! structure-neutral (a molecular rule has no wave vector).
class WaveVectorShift
{
public:
    virtual ~WaveVectorShift() = default;
    //! \f$q\f$ in FRACTIONAL reciprocal coordinates, reduced into [0, 1) per axis -- ONE representative for
    //! every block pair the rule couples, so every pair's \f$\delta\rho(G+q)\f$ is keyed on the same G.
    virtual rvec3_t q() const = 0;
};

} // namespace
