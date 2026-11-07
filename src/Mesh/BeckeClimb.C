// File: BeckeClimb.C  The Becke XC mesh's ONE-ROUND CLIMB: the nested ladder walked to SETTLED at one fixed density.
//
// WHY A CLIMB AT FIXED D (lane B, 2026-11-07; OpenWork "Size every mesh from the seed density").  The geometry pilot judges the Becke mesh
// against the next rung up on each axis and, until now, RE-CONVERGED the SCF after every raise: MnO at eps 1e-5 from the smallest rung is
// 7 raises, so 7 re-convergences.  The seed test (doc/Records/SeedTestReport.md) shows a rung picked at a decent D is the rung the
// converged D picks (5 of 5 systems), and its addendum replays this very rule on banked matrices at one fixed D: 6 of 6 climbs stop on a
// rung that meets eps (4 at the fixed reference's pick, 2 one rung finer).  So the pilot climbs here at the D it has, re-converges ONCE,
// and confirms with one more audit (this function again, with maxRaises 1).  Settings change only at an SCF boundary (pin 31 catch 2).
//
// THE THREE PARTS, kept apart (CLAUDE.md, rigorous vs empirical):
//  - RIGOROUS: none.  The XC mesh has no error bound; every number here is MEASURED (a rung against the rung above, at one D).
//  - THE MODEL: the geometric tail, \c BeckeTailEstimate.  It assumes the rung differences fall geometrically, so a rung's own error is the
//    whole tail d/(1 - r).  It is a shape, not a calibrated constant; where the shape fails (r >= 1) it falls back to d and the next rung
//    measures again.
//  - THE STOPPING RULE: both axes' estimates <= target at the same rung (\c ClimbBeckeLadder).  Never one axis alone: an axis's difference
//    cannot see the other axis's error (LiF nR 30: angular d 1.1e-6 against a true 1.8e-3, which the radial axis's d 1.4e-3 catches).
//  The calibrated part (where the climb STARTS, \c BeckeEntryRung's \c BeckeLadderCalibration) lives in qchem.Mesh.XCPolicy, and changes
//  only how many rungs the climb walks, never where it stops.
module;
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.Mesh.BeckeClimb;
export import qchem.Mesh.XCPolicy;   // MeshParams, BeckeLadderCalibration, BeckeAngularRungAbove / BeckeRadialRungAbove

export namespace qchem::qcMesh
{

//! \brief THE MODEL: a rung's own error from its difference to the rung above, \a d, and the axis's previous difference \a dBelow.  If the
//! differences fall geometrically with ratio r = d/dBelow < 1, the rung's error is the whole tail d/(1 - r).  No \a dBelow (the axis's first
//! measurement), or r >= 1 (not converging), gives d itself: the next rung measures again.  TWO KNOWN PROPERTIES (unit-tested, kept as the
//! pilot has them; doc/Records/SeedTestReport.md addendum):
//!  - without a ratio, d UNDER-reads a geometric tail by (1 - r) -- bounded by 1/(1 - r), and 6 of 6 banked climbs still met eps;
//!  - \a dBelow may come from a rung where the OTHER axis differs, so an unchanged difference reads r ~ 1 and the estimate explodes: an extra
//!    raise (a cost), never an early stop (LiF nR 60 GL 47: d 2.9e-5, estimate 5.0e-4).
//! A same-axis-only ratio with the calibrated factor as fallback was tried and FAILS on NaF (OpenWork, the seed row).
inline double BeckeTailEstimate(double d, std::optional<double> dBelow)
{
    if (dBelow && d<*dBelow) return d/(1.0-d/(*dBelow));
    return d;
}

//! \brief max|H_xc(a) - H_xc(b)| at ONE fixed density, for two meshes, in the normalised basis (pin 31's matrix-element measure).  The
//! caller samples and CACHES: each climb step's reference becomes the next step's rung, so a mesh is sampled once.
using BeckeRungDistance = std::function<double(const MeshParams& a, const MeshParams& b)>;

//! One audit of the climb: \a rung against \a reference, one rung up on \a axis.
struct BeckeClimbStep
{
    MeshParams  rung;
    std::string axis;            //!< "angular" or "radial"
    MeshParams  reference;       //!< the rung above on that axis
    double      delivered=0.0;   //!< max|ΔH_xc| between the two (measured)
    double      tailEstimate=0.0;//!< \c BeckeTailEstimate of it (the model)
    bool        over=false;      //!< the estimate is above the target
};

//! Where and why the climb stopped.
struct BeckeClimbResult
{
    enum class End
    {
        Settled,    //!< both axes meet the target at \c mesh
        AxisEnd,    //!< nothing measured is over, but an axis has no rung above \c mesh, so it is unaudited (\c why names it): not settled, said so
        RaiseCap    //!< \c maxRaises reached: \c mesh is the last raise, NOT audited at this D (a confirming round's single raise)
    };
    MeshParams                  mesh;      //!< the rung the climb stopped at
    std::vector<BeckeClimbStep> steps;     //!< every audit, in order
    int                         raises=0;  //!< rungs climbed
    End                         end=End::Settled;
    std::string                 why;       //!< AxisEnd: the axis that ran out
};

//! \brief THE STOPPING RULE: climb from \a start until both axes' tail estimates meet \a target at one rung.  Each step audits the rung
//! against the rung above on EACH axis (\c BeckeAngularRungAbove, \c BeckeRadialRungAbove), and raises the axis further over.  An axis's
//! estimate uses its own previous difference only when that was measured at a DIFFERENT rung (a re-measure of the same rung has no ratio:
//! d/d ~ 1 made the estimate blow up, NaF geometry 2026-11-02).  This is the pilot's rule (\c SolidCalculation::CalibrateRasters) as it was,
//! lifted out so it can climb at one D and be tested on its own.
//! \param maxRaises  0 = no cap (the ladder's ends bound it: 7 degrees + 5 counts above the bottom).  1 = the pilot's CONFIRMING audit after
//!                   its one re-convergence: measure, and if over, raise once and stop (\c RaiseCap).  THROWS when negative.
BeckeClimbResult ClimbBeckeLadder(const MeshParams& start, double target, const BeckeLadderCalibration& cal, bool siteAdapted,
                                  const BeckeRungDistance& distance, int maxRaises);

} // namespace qchem::qcMesh

namespace qchem::qcMesh
{

static bool SameRung(const MeshParams& a, const MeshParams& b)
{
    return a.nRadial==b.nRadial && a.angularDegree==b.angularDegree && a.angular==b.angular;
}

BeckeClimbResult ClimbBeckeLadder(const MeshParams& start, double target, const BeckeLadderCalibration& cal, bool siteAdapted,
                                  const BeckeRungDistance& distance, int maxRaises)
{
    if (maxRaises<0) throw std::invalid_argument("ClimbBeckeLadder: maxRaises must be >= 0 (0 = no cap)");
    if (!(target>0.0)) throw std::invalid_argument("ClimbBeckeLadder: the target must be > 0");
    BeckeClimbResult out;
    MeshParams cur=start;
    struct Last {MeshParams rung; double d;};
    std::optional<Last> before[2];   // per axis: (rung, difference) at its last measurement
    static const char* const name[2]={"angular", "radial"};
    for (;;)
    {
        const std::optional<MeshParams> up[2]={BeckeAngularRungAbove(cur, cal, siteAdapted), BeckeRadialRungAbove(cur)};
        int worstAxis=-1; double worst=0.0;
        std::string ended;   // an axis with no rung above: unaudited, so it blocks Settled -- but the OTHER axis still climbs (as the pilot did)
        for (int a=0;a<2;++a)
        {
            if (!up[a]) { ended=std::string("the ladder's ")+name[a]+" axis has no rung above "+std::to_string(cur.nRadial)+"/"
                                +std::to_string(cur.angularDegree)+", so that axis is unaudited"; continue; }
            const double d=distance(cur, *up[a]);
            std::optional<double> dBelow;
            if (before[a] && !SameRung(before[a]->rung, cur)) dBelow=before[a]->d;
            const double est=BeckeTailEstimate(d, dBelow);
            before[a]=Last{cur, d};
            out.steps.push_back({cur, name[a], *up[a], d, est, est>target});
            if (est>target && est>worst) { worst=est; worstAxis=a; }
        }
        if (worstAxis<0)
        {
            out.mesh=cur;
            out.end = ended.empty() ? BeckeClimbResult::End::Settled : BeckeClimbResult::End::AxisEnd;
            out.why = ended;
            return out;
        }
        cur=*up[worstAxis]; ++out.raises;
        if (maxRaises>0 && out.raises>=maxRaises) { out.mesh=cur; out.end=BeckeClimbResult::End::RaiseCap; return out; }
    }
}

} // namespace qchem::qcMesh
