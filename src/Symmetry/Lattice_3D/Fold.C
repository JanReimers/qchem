// File: Symmetry/Lattice_3D/Fold.C  Shared orbit-fold primitive for symmetry reductions.
//
// doc/SymmetryUpgradePlan.md §2: ONE fold algorithm serves every reduction surface --
// {k} (BZReduction is re-expressed on FoldGrid), {G} sums (T1), the {r}/Becke quadrature
// (T2), and eventually the GPW stream offsets (T3).  The fold is pure combinatorics: it
// partitions a point set into symmetry orbits (stars) and records, per member, WHICH op
// maps the representative onto it -- the edge datum consumers need for the glide phase
// e^{+2pi i (Um).tau} (G-space), the monomial +/- sign (streams), and diagnostics.
//
// Deliberately NOT here (per the 2026-08-01 review):
//  - No typed containers: each client keeps its native currency ({G}=vector<ivec3_t>,
//    {r}=qcMesh::Mesh, k=IBZMesh) and rebuilds it from the partition.
//  - No weights and no phases: quadrature weights stay real and stay with the client's
//    container; phases are EDGE data recovered from (member, opIndex).
module;
#include <vector>
#include <utility>
export module qchem.Symmetry.Lattice_3D.Fold;
export import qchem.Types;       // ivec3_t, rvec3_t
export import qchem.Matrix3D;    // Matrix3D

export namespace qchem::Symmetry::Lattice_3D
{

//! \brief Spin action carried by a magnetic (Shubnikov) group operation: a plain spatial op
//! leaves the spin channels alone; an op paired with time reversal swaps them (and, on the
//! G-side, conjugates: \f$\tilde\rho(-G)=\tilde\rho(G)^*\f$).  Declared HERE, next to the op
//! type, so every fold consumer is drawn against the spin-native shape from the start
//! (doc/SymmetryUpgradePlan.md §4 tier 4a) -- unpol is the \f$\zeta=0\f$ collapse, not the
//! design case.  Purely spatial folds ignore it: \f$\sigma\f$ acts on FIELD channels, never
//! on the point geometry.
//!
//! NON-COLLINEAR heads-up (plan §4/§9): this enum is itself the COLLINEAR collapse of the
//! general spin action.  A spin-space-group op carries a spin ROTATION -- \f$R_s\in SO(3)\f$
//! on the magnetization 3-vector \f$m(r)\f$, equivalently an SU(2) action on the \f$2\times2\f$
//! spinor density -- with time reversal = \f$m\to-m\f$ + conjugation; {None, Flip} is that
//! action restricted to a fixed quantization axis (spirals/canted orderings need the full
//! rotation).  When non-collinear lands, \f$\sigma\f$ grows into that rotation type; consumers
//! must not bake "two scalar channels" into interfaces beyond the collinear tier.
enum class SpinAction { None, Flip };

//! \brief A generic symmetry operation \f${W|\tau}\f$ for the fold primitives.  \a W is the
//! linear part in the acting space's fractional basis (an integer-valued matrix: the
//! direct-space \f$W\f$ for {r} folds, the reciprocal \f$U\f$ for {G}/{k} folds); \a tau is
//! the fractional translation.  \a tau ACTS on direct-space points (\f$r \to W r + \tau\f$,
//! so a bare matrix cannot fold a non-symmorphic raster), while for {G}/{k} folds it enters
//! only as the edge PHASE \f$e^{+2\pi i\,m\cdot\tau}\f$ recovered via \c Fold::members.
struct SymOp
{
    Matrix3D<double> W;                        //!< Linear part (integer-valued) in the acting space.
    rvec3_t          tau;                      //!< Fractional translation.
    SpinAction       sigma = SpinAction::None; //!< Shubnikov spin action (\f$\sigma\f$), §4 tier 4a.
};

//! \brief The product of a fold: a partition of the raw point set into symmetry orbits,
//! pure combinatorics -- no typed container, no weights, no phases.  "Raw index" is the
//! caller's index into the point set it folded (for \c FoldGrid: the KMesh linear order,
//! ix outer / iz inner).
struct Fold
{
    std::vector<int> owner;      //!< Raw index \f$\to\f$ representative slot.
    std::vector<int> repRaw;     //!< Rep slot \f$\to\f$ raw index of the representative (lowest raw index of its star).
    std::vector<int> starSize;   //!< Rep slot \f$\to\f$ orbit multiplicity \f$n_\mathrm{star}\f$.
    //! Rep slot \f$\to\f$ [(member raw index, op index)]: WHICH op maps rep \f$\to\f$ member --
    //! the hook for the glide phase (G) or monomial sign (streams).  Every orbit member appears
    //! exactly once (the rep included, via its identity/stabilizer op).  An op index of -1 means
    //! no single supplied op realizes the edge (possible only when the op set is not closed under
    //! composition, e.g. a hand-made subset); phase consumers must treat -1 as an error.
    std::vector<std::vector<std::pair<int,int>>> members;

    size_t Reps() const {return repRaw.size();}   //!< Number of orbits (stars).
};

//! \brief Fold arbitrary points under \f$r \to W r + \tau\f$, matching images against the
//! point list within Euclidean tolerance \a tol (the {r}/Becke path).  An op whose image is
//! not in the set is skipped for that point -- folding is merely reduced, never wrong.  No
//! periodic wrap is applied: a torus-periodic point set is the caller's concern (the uniform
//! raster case is exact and belongs to a grid fold).
Fold FoldPoints(const std::vector<rvec3_t>& pts,
                const std::vector<SymOp>& ops, double tol);

//! \brief Fold the periodic grid \f$\{(i+\mathrm{shift})/N\}\f$ under \f$i' = W(i+s) - s
//! \pmod N\f$ -- the EXACT integer path for {k} and {G} grids (raw index = KMesh linear
//! order).  \a tau is NOT applied to the points (on the reciprocal side it is phase-only,
//! recovered from the edge op); an op that does not map the grid onto itself is skipped.
Fold FoldGrid  (const ivec3_t& N, const rvec3_t& shift,
                const std::vector<SymOp>& ops);

//! \brief Fold an explicit G-index list under \f$m \to W m\f$ (exact integer arithmetic;
//! \a W here is the reciprocal G-index map).  \f$|Wm|=|m|\f$, so the Coulomb kernel is
//! constant on each star and a totally symmetric G-sum evaluates at representatives with
//! weight \f$n_\mathrm{star}\f$ (T1).  Images outside the list mean that op is skipped for
//! that point (a symmetric G-ball is closed, so this only trims hand-made sets).
Fold FoldGVectors(const std::vector<ivec3_t>& m,
                  const std::vector<SymOp>& ops);

} // namespace
