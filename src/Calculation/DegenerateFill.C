// File: Calculation/DegenerateFill.C  The facades' policy for a kT = 0 fill of a partly filled DEGENERATE level that nothing stated
// (D-DEGEN-FILL, user rulings 2026-11-08), and its one report row.
//
// The fill RULE is an injected object (qchem.ElectronConfiguration.OccupationPolicy's occupancy axis): IntegerOccupancy fills
// such a level orbital by orbital, SharedDegenerateOccupancy gives every member the same share (the kT -> 0 Fermi limit: D
// invariant to the eigensolver's rotation inside the level).  The caller states it (SCFParams::DegenerateLevels, deck key
// degenerateLevelFill); when it states Auto, the FACADE chooses here -- the one place the choice lives.  Either way the
// calculation files ONE `choose` row naming the rule that fills and why (user: no further warning rows, "C will confuse users").
module;
#include <optional>
#include <string>
export module qchem.DegenerateFillReport;
export import qchem.SCFParams;
import qchem.Policy;
import qchem.Reporting;

export namespace qchem
{

//! \brief THE DEGENERATE-LEVEL RULE of a kT = 0 run: the stated one, or (on \c Auto) the facade's \a autoChoice; files ONE `choose`
//! row naming the rule and why (\a report: once per calculation).  \a why = the facade's reason for its default.
inline SCFParams ResolveDegenerateLevelFill(const SCFParams& asked, SCFParams::DegenerateLevelFill autoChoice, const std::string& why,
                                            bool report)
{
    SCFParams p=asked;
    const bool stated = p.DegenerateLevels!=SCFParams::DegenerateLevelFill::Auto;
    if (!stated) p.DegenerateLevels=autoChoice;
    if (report && p.SmearingkT<=0.0)
        qchem::policy::Report("DegenerateLevelFill",
            p.DegenerateLevels==SCFParams::DegenerateLevelFill::Shared
                ? "SharedDegenerateOccupancy: a partly filled degenerate level is shared equally among its orbitals (degenerateLevelFill: Shared)"
                : "IntegerOccupancy: a partly filled degenerate level is filled orbital by orbital (degenerateLevelFill: OrbitalByOrbital)",
            qchem::policy::Verdict::Choose, stated ? std::string("stated by the caller (degenerateLevelFill)") : why);
    return p;
}

} // namespace qchem
