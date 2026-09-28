// File: Response/Reference.C  The UNPERTURBED STATE a linear response linearises about, and its
// independent-particle response R0 (doc/LinearResponsePlan.md §1–§3, stage R0).
//
// WHAT IT ANSWERS.  Given a first-order Fock change on the block pairs a wave-vector shift q couples,
// (k, σ) -> (k+q, σ), what first-order density matrix does a system of INDEPENDENT particles produce?
//     δD_{mn} = W_{mn} δF_{mn},   W_{mn} = g (f_{nk} - f_{m,k+q}) / (ε_{nk} - ε_{m,k+q})
// in the orbital (MO) basis -- a dense SUM OVER STATES.  A Gaussian basis diagonalises the whole space at
// every k, so every virtual is at hand: no Sternheimer solve, no P_c projector, no k+q NSCF run (the plane-
// wave machinery of QE's `sternheimer_kernel` exists because a PW code cannot afford the empty states).
// A future large-cutoff PW code joins as another concrete of R0, not as a change here (§4b).
//
// THREE CONTRACTS, each a ruling of the plan's review round:
//  * THE WEIGHT COMES FROM THE RUN'S OWN OCCUPANCY RULE (E1).  The reference is built over the SAME rule the
//    ground state filled with (qchem::MakeOccupancyRule from the run's OccupationConfig), so integer vs
//    Fermi is decided once, for the ground state and its response alike.
//  * THE RESPONSE GAP IS MEASURED, AND GATED WHERE IT MUST BE (E1, review comment 2).  An integer rule's
//    weight diverges as a coupled occupied/empty pair closes.  Crystal_EC's insulator mode fills a FIXED
//    count per k-block with no aufbau across blocks, so across a (k, k+q) pair nothing orders the levels:
//    an occupied band at k+q can sit ABOVE an empty one at k.  Gap() returns the smallest
//    Δ = ε_empty - ε_occ over the pairs q couples and FAILS if it is <= 0 (inverted) or not resolved above
//    the reference's own eigenvalue noise.  A small POSITIVE gap is physics (a large χ0), never flagged.
//  * THE PARTNER BLOCK MUST BE STORED (D5): a full, unreduced k-mesh.  An IBZ-reduced reference (star > 1)
//    is refused at construction -- a broken precondition, so it throws.
module;
#include <limits>
#include <memory>
#include <string>
#include <vector>
export module qchem.Response.Reference;
export import qchem.Symmetry.Irrep;
export import qchem.Symmetry.Lattice_3D.BlochQN;               // MeshShift (the lattice selection rule, and the q-mesh)
export import qchem.Symmetry.SelectionRule;                     // which block a perturbation couples to which (S1)
export import qchem.ElectronConfiguration.OccupationPolicy;    // OccupancyRule
export import qchem.Outcome;
export import qchem.Types;

export namespace qchem::Response
{
using Symmetry::Lattice_3D::MeshShift;
using Symmetry::SelectionRule;
using cmat_t=mat_t<dcmplx>;

//! \brief One stored block of the unperturbed state: a (k, σ) Bloch block of a FULL mesh.  A construction
//! VALUE (like OccupationConfig): the wave-function adapter fills it, and so can a unit test's model.
struct ReferenceBlock
{
    Irrep  irrep;          //!< the block's identity: (k, σ); \c irrep.sym is a BlochQN of star 1
    double w=0.0;          //!< BZ weight (\f$1/N_{\rm mesh}\f$ on a full mesh)
    double g=1.0;          //!< level capacity: the spin degeneracy (2 for a folded doublet, 1 per channel)
    rvec_t e;              //!< EVERY orbital's eigenvalue (Hartree), occupied and virtual, in stored order
    rvec_t f;              //!< the matching FRACTIONAL occupations, occupation/g, in [0,1]
    int    reservoir=0;    //!< blocks sharing one chemical potential share an id -- the q=0 Fermi shift δμ
};

//! Matrices on the block pairs ONE selection rule couples, in the orbital basis: \c m[b] is
//! \f$n_{\rm orb}({\rm bra})\times n_{\rm orb}(b)\f$ for ket block \a b and its partner bra block (k+q for a
//! lattice shift, \a b itself for an \c Invariant perturbation).
struct BlockPairs
{
    std::vector<cmat_t> m;
};

//! WHY a response could not be computed.  A VALUE (an Outcome's error), never a throw: the caller chose
//! the state and the mesh and can choose again.
struct ResponseFailure
{
    enum class Why {Inverted, Unresolved, Incommensurate};
    Why         why=Why::Unresolved;
    std::string detail;
};

//! The measured response gap of one q, and what it implies.
struct ResponseGap
{
    double gap  =std::numeric_limits<double>::infinity();   //!< min Δ over the gated coupled pairs (Hartree); inf when none is gated
    double noise=0.0;                                         //!< the reference's eigenvalue noise δε (Hartree); NaN = unmeasured
    //! The relative bound on χ0 it implies, \f$\delta\chi_0/\chi_0\lesssim\delta\varepsilon/\Delta_{\min}\f$.
    double Bound() const {return noise/gap;}
};

class Reference
{
public:
    //! \a rule: the run's own occupancy rule (\c MakeOccupancyRule of the run's configuration).
    //! \a eigenNoise: how well the reference's eigenvalues are converged (Hartree) -- MEASURED by the
    //! caller (the SCF's final [F,D] commutator), never a tunable: it is the threshold below which a gap
    //! is not resolved.  NaN = UNMEASURED (a recipe that computes no commutator): the gate then checks the
    //! gap's SIGN only, and every report says "unmeasured" rather than printing a false 0.
    //! THROWS if a block is not a star-1 Bloch block (D5).
    Reference(std::vector<ReferenceBlock> blocks, std::unique_ptr<OccupancyRule> rule, double eigenNoise);

    size_t       NumBlocks  ()         const {return itsBlocks.size();}
    size_t       NumOrbitals(size_t b) const {return itsBlocks[b].e.size();}
    double       Weight     (size_t b) const {return itsBlocks[b].w;}
    const Irrep& BlockIrrep (size_t b) const {return itsBlocks[b].irrep;}
    double       EigenNoise ()         const {return itsNoise;}

    //! The q-mesh of \a Nq divisions on this reference's k-mesh -- FAILS if incommensurate.
    Outcome<std::vector<MeshShift>,ResponseFailure> QMesh(ivec3_t Nq) const;
    //! For every ket block, the index of the ONE bra block \a rule couples it to (same spin): k+q for a
    //! \c MeshShift, itself for an \c Invariant perturbation.  THROWS unless exactly one partner is stored:
    //! for a lattice, the full-mesh precondition (D5) was broken, which no caller can repair here.
    //! (Exactly one: a point-group product rule that couples one ket to SEVERAL bra irreps generalises the
    //! BlockPairs layout, when it is built.)
    std::vector<size_t> Partners(const SelectionRule& rule) const;
    //! The RESPONSE GAP over the pairs \a rule couples (see the file header): FAILS on an inverted or
    //! unresolved pair when the rule needs a resolved gap; otherwise the measured gap (inf for a smeared rule).
    Outcome<ResponseGap,ResponseFailure> Gap(const SelectionRule& rule) const;
    //! \brief R0: a first-order Fock change \a dF on \a rule's block pairs (orbital basis) -> the independent-
    //! particle first-order density matrix.  When every block is its OWN partner (q = 0, or a totally
    //! symmetric perturbation) it adds the Fermi-level shift of each RESERVOIR,
    //! \f$\delta\mu_r=\sum_{b\in r}w_b\sum_n D_n\delta F_{nn}/\sum_{b\in r}w_b\sum_n D_n\f$ with \f$D_n=g f'_n\f$,
    //! so every reservoir keeps its electron count (\f$\delta f_n=f'_n(\delta\varepsilon_n-\delta\mu)\f$).
    //! An integer rule has \f$f'\equiv0\f$, so it shifts nothing.
    //! \warning Call Gap(rule) first: an ungated integer pair throws inside the weight (a broken invariant).
    BlockPairs ApplyR0(const SelectionRule& rule, const BlockPairs& dF) const;
    //! \f$\sum_b w_b\sum_{mn}\bar a_{b,mn}\,x_{b,mn}\f$ -- the BZ-weighted pairing of an operator with a
    //! first-order density: an expectation value's first-order change, per unit cell.
    dcmplx Contract(const BlockPairs& a, const BlockPairs& x) const;

private:
    std::vector<ReferenceBlock>    itsBlocks;
    std::unique_ptr<OccupancyRule> itsRule;
    double                         itsNoise;
};

} // namespace
