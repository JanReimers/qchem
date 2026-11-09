// File: Symmetry/StarGroup.C  The group a run star-averages its density under, and the face a field offers it.
//
// ★ THE OPS STAY WITH THE GROUP (V-SYM-SEQUESTER-OPS, first piece; user 2026-11-08/09).  A k-composite density used to hold a
// vector of reciprocal ops and hand it to whoever averaged (the projector for ρ̃, the fit basis -- which stored its OWN copy --
// for the raster ρ).  Now the composite holds the GROUP and asks it: `group.StarAverage(field)`.  The group iterates its own
// ops; the FIELD knows how it transforms under one op (a G-map permutes its indices with the glide phase, a raster permutes its
// voxels after an exact FFT shift).  So what the basis KNOWS (how its functions transform) is split from what the run CHOOSES
// (which group is imposed), and no op list travels.
//
// ★ STRUCTURE-NEUTRAL (user 2026-11-09: "for the purpose of interface sanity we can think of molecules as k=0 systems").  The
// ops are \c Symmetry::SymOp \f$\{W|\tau\}\f$ in the acting space's own basis: fractional for a crystal, Cartesian with
// \f$\tau=0\f$ for a molecule's point group.  Only a crystal needs a star-average today (an IBZ k-sum builds a partial density;
// a molecule's SALC-blocked D is invariant by construction), but a molecule is the k=0 case of the same interface, not a
// different one.
//
// Not the tier-1 \c Representation<T> hub: that one acts on MATRICES (\f$M^\dagger OM\f$, \f$MXM^\dagger\f$); a star-average
// acts on a one-index FIELD, which can reach outside any basis's function list (a ρ̃ map keyed by G, a raster), so it is a face
// of its own (user ruling B, 2026-11-09).
module;
#include <memory>
#include <vector>
export module qchem.Symmetry.StarGroup;
export import qchem.Symmetry.SymOp;

export namespace qchem::Symmetry
{

//! \brief A field the group can star-average: it adds its own IMAGE under one op into an accumulator it keeps.
//!
//! The protocol is the group's (\c StarGroup::StarAverage): \c StartStar once, \c AddImage once per op (in the group's op
//! order), then \c FinishStar with the weight \f$1/|G|\f$, which replaces the field by the weighted image sum.  HOW the sum is
//! formed is the field's: the raster accumulates each image as it comes and scales once; the G-map DEFERS, collecting the
//! ops and sweeping its own entries once at \c FinishStar.  Each keeps the arithmetic order its op-vector average had, so
//! moving the ops into the group changed no bits.
class StarField
{
public:
    virtual ~StarField() = default;
    virtual void StartStar ()=0;                      //!< zero the accumulator
    virtual void AddImage  (const SymOp& g)=0;        //!< the field's image under \a g joins the sum
    virtual void FinishStar(double weight)=0;         //!< field <- \a weight * accumulator
};

//! \brief A finite group of \c SymOp s that star-averages fields: \f$f\leftarrow|G|^{-1}\sum_{g\in G}g\cdot f\f$.
//!
//! A PROJECTOR only over a closed group -- the imposition gate checks closure before one is built (ResolveImposition).  The
//! EMPTY group is the trivial \f$\{E\}\f$: its average is an exact no-op, so no caller asks whether symmetry was imposed.
class StarGroup
{
public:
    StarGroup() = default;                                            //!< the trivial group {E}
    explicit StarGroup(std::vector<SymOp> ops) : itsOps(std::move(ops)) {}

    bool   IsTrivial() const {return itsOps.empty();}
    size_t size     () const {return itsOps.size();}                  //!< |G| (0 for the trivial group, which stores no E)

    //! Star-average \a f over my ops, in place (a no-op for the trivial group).
    void StarAverage(StarField& f) const
    {
        if (itsOps.empty()) return;
        f.StartStar();
        for (const SymOp& g : itsOps) f.AddImage(g);
        f.FinishStar(1.0/double(itsOps.size()));
    }

private:
    std::vector<SymOp> itsOps;   //!< the group's ops, in the acting space's own basis (\c SymOp)
};

} // namespace qchem::Symmetry
