// File: Calculation/DegenerateFill.C  The facades' policy for a kT = 0 fill of a partly filled DEGENERATE level that nothing stated
// (D-DEGEN-FILL, user rulings 2026-11-08), and its one report row.
//
// The fill RULE is an injected object (qchem.ElectronConfiguration.OccupationPolicy's occupancy axis): IntegerOccupancy fills
// such a level orbital by orbital, SharedDegenerateOccupancy gives every member the same share (the kT -> 0 Fermi limit: D
// invariant to the eigensolver's rotation inside the level).  The caller states it (SCFParams::DegenerateLevels, deck key
// degenerateLevelFill); when it states Auto, the FACADE chooses here -- the one place the choice lives.  Either way the
// calculation files ONE `choose` row naming the rule that fills and why (user: no further warning rows, "C will confuse users").
module;
#include <algorithm>
#include <optional>
#include <sstream>
#include <string>
export module qchem.DegenerateFillReport;
export import qchem.SCFParams;
import qchem.Policy;
import qchem.Reporting;

export namespace qchem
{

//! \brief What a declined direct minimiser falls back to, in words, for a policy row: the mixer \a resolved names (Auto already
//! resolved by the iterator, or \c Linear on a molecule, whose mixer is linear D-mixing) and the relax the iterator will run it at --
//! a linear (adaptive) mixer under the \c DeclinedDirectMinRelax ceiling (CLAUDE.md: 0.2-0.3, never 1.0), every history mixer
//! (Pulay on D or on rho~) and Kerker at its own.  Pulay-D is deliberately NOT capped (lane A, 2026-11-09, Shared + GDM declining
//! at every step, stated Pulay-D depth 8 at alpha 0.45/0.3/0.25/0.2: Si triplet box 12/11/11/11 iterations, MnO dimer 67 / cap /
//! cap / 57 -- the step is not what converges a Pulay history, and a 0.25 ceiling would have failed the MnO case).
inline std::string DescribeDeclinedFallback(const SCFParams& resolved)
{
    using M=SCFParams::Mixing;
    const double cap=std::min(resolved.StartingRelaxRo, resolved.DeclinedDirectMinRelax);
    std::ostringstream os;
    if (resolved.Mixer==M::Pulay && resolved.KerkerG0<=0.0)
        os<<"Pulay D-mixing (depth "<<resolved.PulayDepth<<") at relax "<<resolved.StartingRelaxRo;
    else if (resolved.Mixer==M::Pulay)  os<<"Kerker + Pulay density mixing at relax "<<resolved.StartingRelaxRo;
    else if (resolved.Mixer==M::Kerker) os<<"Kerker density mixing at relax "<<resolved.StartingRelaxRo;
    else                                os<<"linear D-mixing at relax "<<cap;
    return os.str();
}

//! What the facade knows about the run's accelerator that bears on the fill: a direct minimiser (GDM) rotates an INTEGER-occupation
//! determinant, so it declines a density whose degenerate level is shared.
struct DegenerateFillContext
{
    bool        directMin = false;   //!< the run's accelerator is GDM (stated, or resolved to it)
    std::string fallback;            //!< \c DescribeDeclinedFallback of the run's mixer: what runs when GDM declines
};

//! \brief THE DEGENERATE-LEVEL RULE of a kT = 0 run: the stated one, or (on \c Auto) the facade's \a autoChoice; files ONE `choose`
//! row naming the rule and why (\a report: once per calculation).  \a why = the facade's reason for its default.
//!
//! GDM (user 2026-11-08: "Many users will have no idea that GDM only works for integers. So we need to guide them"): on \c Auto a
//! GDM run takes the integer fill, and the row says why and how to get the symmetric state; a GDM run that STATES \c Shared gets a
//! `warn` row saying what will happen instead (GDM declines a shared level, the named fallback runs), so a "GDM run" that never ran
//! GDM is no surprise.
inline SCFParams ResolveDegenerateLevelFill(const SCFParams& asked, SCFParams::DegenerateLevelFill autoChoice, const std::string& why,
                                            bool report, const DegenerateFillContext& ctx = {})
{
    using F=SCFParams::DegenerateLevelFill;
    SCFParams p=asked;
    const bool stated = p.DegenerateLevels!=F::Auto;
    std::string reason = stated ? std::string("stated by the caller (degenerateLevelFill)") : why;
    if (!stated && ctx.directMin)
    {
        autoChoice=F::OrbitalByOrbital;
        reason="GDM works only on integer occupations, so a partly filled degenerate level is filled orbital by orbital (a "
               "symmetry-broken state); state degenerateLevelFill: Shared for the symmetric state, and GDM will then hand over to "
               "the fallback mixer";
    }
    if (!stated) p.DegenerateLevels=autoChoice;
    if (!report || p.SmearingkT>0.0) return p;
    qchem::policy::Report("DegenerateLevelFill",
        p.DegenerateLevels==F::Shared
            ? "SharedDegenerateOccupancy: a partly filled degenerate level is shared equally among its orbitals (degenerateLevelFill: Shared)"
            : "IntegerOccupancy: a partly filled degenerate level is filled orbital by orbital (degenerateLevelFill: OrbitalByOrbital)",
        qchem::policy::Verdict::Choose, reason);
    if (ctx.directMin && p.DegenerateLevels==F::Shared)
        qchem::policy::Report("DegenerateLevelFillUnderGDM",
            "GDM will DECLINE a shared degenerate level and the run takes "+ctx.fallback+" instead",
            qchem::policy::Verdict::Warn,
            "GDM rotates an integer-occupation determinant, and a partly filled degenerate level shared equally is fractional: if the "
            "run meets one, GDM declines every step and the run is mixed fixed-point steps, not GDM.  State degenerateLevelFill: "
            "OrbitalByOrbital for a GDM run (a symmetry-broken state)");
    return p;
}

} // namespace qchem
