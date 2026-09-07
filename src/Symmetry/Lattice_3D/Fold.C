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
#include <cassert>   // the SymmetrizeValuesSigned edge-op guard
export module qchem.Symmetry.Lattice_3D.Fold;
export import qchem.Types;       // ivec3_t, rvec3_t
export import qchem.Matrix3D;    // Matrix3D
export import qchem.Symmetry.SymOp;   // SymOp / SpinAction (structure-neutral; aliased below)

export namespace qchem::Symmetry::Lattice_3D
{

// SpinAction and SymOp were promoted to the structure-NEUTRAL qchem.Symmetry.SymOp (R2.17.3): the
// operation {W|tau} is not a crystal concept -- a molecule's point group is the tau=0 case -- and living
// in this folder made every signature that took one look crystal-specific.  Aliased here so that every
// existing Symmetry::Lattice_3D::SymOp / ::SpinAction spelling keeps working unchanged; the fold
// primitives below are, and remain, lattice-specific (they fold on the fractional torus).
using qchem::Symmetry::SpinAction;
using qchem::Symmetry::SymOp;

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
//! periodic wrap is applied: a torus-periodic point set takes \c FoldPointsPeriodic.
Fold FoldPoints(const std::vector<rvec3_t>& pts,
                const std::vector<SymOp>& ops, double tol);

//! \brief Is \a ops CLOSED under composition -- i.e. actually a GROUP?
//!
//! ★ WHY THIS EXISTS (2026-09-07, doc/SymmetryUpgradePlan.md "SUPERCELLS").  Every star-average in the
//! tree assumes its op set is a group, because that is what makes the average a PROJECTOR
//! (\f$P^2=P\f$).  Nothing checked it, and a set that is not closed corrupts a density silently.
//! MEASURED: `SpaceGroup::Detect` on a Si 2x2x2 SUPERCELL returns 48 ops of which **864 of the 2304
//! products fall outside the set** -- because Detect assumes a primitive cell (one \f$\tau\f$ coset per
//! \f$W\f$) and on an 8-fold non-primitive cell \f$\tau\f$ is only defined modulo the internal
//! translations, so arbitrary per-\f$W\f$ representatives do not compose.  The run diverged from its
//! first Fock and three SCF explanations were tested before anyone suspected the group.
//!
//! ⚠ **CLOSURE, NOT COMPLETENESS, IS THE INVARIANT.**  A proper SUBGROUP is perfectly legal here -- the
//! imposed path deliberately hands the legacy faces only the \f$\sigma=\f$None subgroup, and a partial
//! projector merely folds less.  What is never legal is a set that is not a group at all.
//!
//! Composition is \f$\{W_1|\tau_1\}\{W_2|\tau_2\}=\{W_1W_2\,|\,W_1\tau_2+\tau_1\}\f$ with
//! \f$\tau\f$ compared modulo 1, and the Shubnikov \f$\sigma\f$ multiplying as \f$\pm1\f$
//! (Flip\f$\circ\f$Flip = None) -- a magnetic set must close in \f$\sigma\f$ too.
//! \note \f$O(|G|^3)\f$ naive, on groups of tens of ops, ONCE per run: microseconds.  Cheap enough that
//! there is no argument for not doing it.
bool IsClosedGroup(const std::vector<SymOp>& ops, double tol=1e-6);

//! \brief Throw unless \a ops is closed (see \c IsClosedGroup); \a who names the imposition site.
//!
//! THROWS rather than asserts, and deliberately: every production run is Release, where \c assert is
//! compiled out -- so an assert would have caught exactly none of the runs this defect actually spoiled.
//! Same reasoning as \c RequireSiteBlocks in UnitCell.C: a silently wrong answer is worse than a stopped
//! run.  A broken group is a broken INVARIANT no caller can act on, so it throws rather than returning an
//! Outcome (CLAUDE.md).
void RequireClosedGroup(const std::vector<SymOp>& ops, const char* who, double tol=1e-6);

//! \brief The TORUS variant of \c FoldPoints for cell-periodic point sets (the periodic Becke
//! mesh, T2): points are FRACTIONAL coordinates, images \f$W f + \tau\f$ are matched modulo the
//! lattice (each component reduced into \f$[0,1)\f$, distance = the torus metric).  Input points
//! need not be pre-wrapped.
Fold FoldPointsPeriodic(const std::vector<rvec3_t>& pts,
                        const std::vector<SymOp>& ops, double tol);

//! \brief The INVARIANCE checker (§6 T2 precondition): per op, how many points' images
//! \f$W f + \tau\f$ (torus metric) fail to land back on the point set.  All zeros == the set is
//! invariant under every op -- the precondition for pointwise star-averaging (\c SymmetrizeValues
//! as a projector) and for weight folding to hit full multiplicity.
std::vector<int> CountUnmappedPeriodic(const std::vector<rvec3_t>& pts,
                                       const std::vector<SymOp>& ops, double tol);

//! \brief Pointwise STAR-AVERAGE of per-point values over a fold's orbits: replace every
//! member's value by its orbit mean -- the real-space {r} sibling of the G-space density
//! star-average.  This equals the exact group projector \f$(Pf)(r)=\tfrac1{|G|}\sum_g f(g\,r)\f$
//! ONLY when \a ops form a full group AND the point set is invariant under it (then each
//! point's image multiset covers its orbit uniformly, \f$|G|/|\mathrm{star}|\f$ copies each);
//! on a partial op set it is merely a smoothing, so consumers must fold under a genuine
//! (sub)group -- the §3 policy's job.
template <class V> inline void SymmetrizeValues(const Fold& f, V& vals)
{
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        double s = 0.0;
        for (auto [m, o] : f.members[r]) s += vals[m];
        s /= double(f.members[r].size());
        for (auto [m, o] : f.members[r]) vals[m] = s;
    }
}

//! \brief The SIGNED star-average for a field of ODD spin parity under a MAGNETIC (Shubnikov)
//! group -- S2 of doc/SymmetryUpgradePlan.md §7 step 7.  The channel pair diagonalizes the spin
//! action: the total \f$\rho_\uparrow+\rho_\downarrow\f$ is EVEN under \c Flip (take the plain
//! \c SymmetrizeValues over the same fold), while the magnetization \f$m=\rho_\uparrow-
//! \rho_\downarrow\f$ picks up the character \f$\chi(g)=-1\f$ on \f$\sigma=\f$\c Flip ops:
//! \f$(Pm)(r)=\tfrac1{|M|}\sum_g \chi_g\,m(g\,r)\f$.  \a sigmas is the per-op spin action,
//! PARALLEL to the op list the fold was built under (the edge opIndex indexes it); rep value
//! \f$s=\tfrac1n\sum \chi\,v\f$, member scatter \f$v=\chi\,s\f$ (\f$\chi=\pm1\f$ is its own
//! inverse).
//!
//! PRECONDITIONS (beyond the plain version's full-group + invariant-set): a point whose
//! STABILIZER contains a \c Flip op carries \f$m\equiv0\f$ exactly (the projector annihilates the
//! odd field there), and the single-edge orbit mean CANNOT see that -- callers must zero \a vals
//! at the points \c FlipFixedPointsPeriodic flags BEFORE trusting this as the projector.  (MnO
//! AFM-II: the O sites are fixed by the Mn-sublattice-swapping inversion, and m(O)=0 is exactly
//! the physics.)  An edge with no mapping op (-1) cannot carry a character and asserts.
template <class V> inline void SymmetrizeValuesSigned(const Fold& f,
                                                      const std::vector<SpinAction>& sigmas, V& vals)
{
    auto chi=[&](int o)->double
    {
        assert(o>=0 && "SymmetrizeValuesSigned: an edge with no mapping op cannot carry a character");
        return sigmas[size_t(o)]==SpinAction::Flip ? -1.0 : 1.0;
    };
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        double s = 0.0;
        for (auto [m, o] : f.members[r]) s += chi(o)*vals[m];
        s /= double(f.members[r].size());
        for (auto [m, o] : f.members[r]) vals[m] = chi(o)*s;
    }
}

//! \brief The ODD-field fixed-point AUDIT for \c SymmetrizeValuesSigned: flag every point some
//! \f$\sigma=\f$\c Flip op maps onto ITSELF (torus metric).  The exact projector annihilates an
//! odd field at such points, so the caller zeroes them; geometry-fixed, so run it once per mesh.
std::vector<char> FlipFixedPointsPeriodic(const std::vector<rvec3_t>& pts,
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
