// File: Mesh/FoldedMesh.C  A quadrature mesh TOGETHER WITH the orbit partition it is invariant under.
//
// ★ WHY THIS IS ONE TYPE AND NOT TWO THINGS PASSED SIDE BY SIDE (user, 2026-09-08).
//
// The pointwise star-average -- replace every value by its orbit mean -- is an EXACT PROJECTOR **only
// because the mesh is invariant under the ops the fold was built from** (doc/SymmetryUpgradePlan.md §6a
// W1).  That is a joint property of a mesh AND a fold, and until now nothing expressed it: the two
// travelled side by side in `BasisSet::FitQuadrature` -- declared in a FITTING header, which is neither
// mesh nor symmetry business -- and the consumer wired one to the other by hand, asserting sizes at the
// point of use.  A caller could hand over a fold built for a different mesh and find out at the first
// out-of-range index, or not at all.
//
// Here the pairing is a CLASS INVARIANT, checked once at construction, and the star-average is a member
// rather than something every consumer re-writes.  It is the same move as keeping an integrator's forward
// and adjoint in one object: when two things are only correct together, make them one thing.
//
// ⚠ THIS IS WHY qcMesh LINKS qcSymmetry.  `Fold` is pure INDEX bookkeeping -- owner / repRaw / starSize /
// members, no geometry -- and `SymmetrizeValues` is a template over the value container, so the dependency
// is on a partition type, not on group theory.  qcSymmetry imports nothing from qcMesh, so this closes no
// cycle (checked both directions, 2026-09-08).
module;
#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.Mesh.Folded;
export import qchem.Mesh;                       // qcMesh::Mesh
export import qchem.Symmetry.Lattice_3D.Fold;   // Fold + SymmetrizeValues / SymmetrizeValuesSigned
export import qchem.Symmetry.Spin;              // Symmetry::SpinAction (the Shubnikov tags)

export namespace qchem::qcMesh
{

//! \brief A quadrature mesh plus the orbit partition it is invariant under -- and, on a MAGNETICALLY
//! imposed run, the Shubnikov spin tags that go with it.
//!
//! DEFAULT-CONSTRUCTED = A FREE RUN: no mesh, no fold, and every symmetrisation below is an exact no-op.
//! That is the common case and it must stay silent -- no consumer should have to ask "was symmetry
//! imposed?" before calling.
class FoldedMesh
{
public:
    using Fold       = Symmetry::Lattice_3D::Fold;
    using SpinAction = Symmetry::SpinAction;

    FoldedMesh() = default;

    //! \brief The checked constructor -- THE reason this type exists.
    //!
    //! \a fold's indices must address \a mesh's points, and \a flipFixed (when present) must be one flag
    //! per point.  A fold built for a different mesh is a composition error: it produces silently wrong
    //! averages where the indices happen to be in range, so it THROWS rather than asserting -- an assert
    //! is compiled out under NDEBUG, which is exactly where production runs live.
    FoldedMesh(std::shared_ptr<const Mesh> mesh, Fold fold,
               std::vector<SpinAction> sigmas={}, std::vector<char> flipFixed={});

    //! The mesh.  Null on a default-constructed (free-run) bundle -- ask before dereferencing.
    const std::shared_ptr<const Mesh>& GetMesh() const {return itsMesh;}
    bool   HasFold  () const {return !itsFold.owner.empty();}   //!< is a partition present at all?
    bool   IsMagnetic() const {return !itsSigmas.empty();}      //!< are Shubnikov spin tags present?
    size_t size     () const {return itsMesh ? itsMesh->size() : 0;}     //!< mesh points
    //! Orbits under the fold -- \c size() when there is no fold, since then every point is its own star.
    //! Reported as the FOLD FACTOR's denominator, so the no-fold answer has to be the honest one.
    size_t NumOrbits() const {return HasFold() ? itsFold.Reps() : size();}
    //! How many ops the Shubnikov tags cover (0 = not magnetically imposed).  Reporting only.
    size_t NumSpinOps() const {return itsSigmas.size();}

    //! \brief STAR-AVERAGE a per-point field in place: the exact orbit projector on an invariant mesh.
    //! No fold => exact no-op, so no caller asks whether symmetry was imposed.  REAL-space, so it
    //! PRESERVES \f$\rho\ge0\f$.
    void Symmetrize(rvec_t& f) const;

    //! \brief The MAGNETIC sibling: project the \f$(\rho,m)\f$ PAIR, which is what diagonalises \f$\sigma\f$
    //! -- \f$\rho\f$ EVEN under the orbit mean, \f$m\f$ ODD under the \f$\chi\f$-signed one, with the
    //! flip-fixed entries of \f$m\f$ zeroed first (Shubnikov S3, doc/SymmetryUpgradePlan.md §7).
    //! No spin tags => grey/free semantics => each channel averaged independently.
    void SymmetrizeSpin(rvec_t& rho, rvec_t& m) const;

    //! \brief The Shubnikov spin actions, parallel to the op list the fold was built under.  READ-only,
    //! and for OBSERVATION (a gate counting flips, a report line) -- the projection itself is
    //! \c SymmetrizeSpin above, which is where the tags and the fold are used together correctly.
    const std::vector<SpinAction>& GetSpinOps() const {return itsSigmas;}
    //! Is mesh point \a g fixed by some \f$\sigma\f$=Flip op?  \f$m\f$ must vanish there exactly.
    bool   IsFlipFixed(size_t g) const {return g<itsFlipFixed.size() && itsFlipFixed[g];}
    size_t NumFlipFixed() const {return itsFlipFixed.size();}   //!< 0 = no odd-field audit was supplied

    //! The raw partition, for the few consumers that must walk orbits themselves (the orbit-consistency
    //! filter).  Prefer the two members above: they are why this type exists.
    const Fold& GetFold() const {return itsFold;}

private:
    std::shared_ptr<const Mesh> itsMesh;
    Fold                        itsFold;      //!< its orbit partition ({} = no star-averaging)
    //! Shubnikov S3: the per-op spin actions PARALLEL to the op list the fold was built under (the edge
    //! opIndex indexes it), and the odd-field zero flags (points some \f$\sigma\f$=Flip op maps onto
    //! themselves, where \f$m\f$ must vanish exactly).  Both empty = grey/free semantics.
    std::vector<SpinAction>     itsSigmas;
    std::vector<char>           itsFlipFixed;
};

} // namespace
