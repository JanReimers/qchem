// File: WaveFunction/Imp/EigenRows.C  The pure row builders behind the level tables (D10).
module;
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>
#include <string>
#include <vector>
module qchem.WaveFunction.EigenRows;

namespace qchem::WaveFunction
{
using Orbitals::EnergyLevel;
using Orbitals::EnergyLevels;

std::string OccupationCell(double occ, int degen)
{
    std::ostringstream os;
    // Integer occ (gapped insulator) unchanged; FRACTIONAL (Fermi-smeared) shown with decimals.  setprecision(0)
    // alone rounded a smeared 0.996 to "1/1" and 0.004 to "0/1" -- an integer configuration for a run whose own
    // trace column was flagging partial occupancy every iteration.
    if (std::fabs(occ-std::round(occ)) < 1e-6) os << std::fixed << std::setprecision(0) << occ;
    else                                       os << std::fixed << std::setprecision(2) << occ;
    os << "/" << degen;
    return os.str();
}

static std::string Fixed8(double x) { std::ostringstream os; os << std::fixed << std::setprecision(8) << x; return os.str(); }

static std::string LabelOf(const EnergyLevel& el) { std::ostringstream os; os << el.qns.n << *el.qns.sym; return os.str(); }

std::vector<EigenRow> UnpolarizedEigenRows(const EnergyLevels& els)
{
    // The HIGHEST occupied energy -- the honest end of the table.  Occupations are monotonic in energy under one
    // mu, so for AUFBAU this is the level before the first empty one.  It is NOT monotonic under MOM: a
    // character-pinned run can leave a level EMPTY well BELOW an occupied one (the hole the 0h guard watches for),
    // and a plain `break` at the first empty level then truncates the table exactly AT the anomaly.  So run to the
    // highest OCCUPIED level, never stopping short of it (measured on MnO: a -1.29 Ha EMPTY level, invisible in a
    // table that happily printed a +0.75 Ha virtual).
    double eHomo=-1e300;
    for (const auto& [e,el] : els) if (el.occ >= 1e-6) eHomo=e;
    std::vector<EigenRow> rows;
    for (const auto& [e,el] : els)
    {
        // Stop past the frontier by OCCUPATION, not energy sign.  The old `e>0.0` cutoff is a MOLECULAR idiom; in a
        // SOLID the energy zero is arbitrary, so the Fermi level -- and every occupied level -- can be POSITIVE.
        if (el.occ < 1e-6 && e > eHomo) break;
        rows.push_back({OccupationCell(el.occ, el.degen), Fixed8(e), LabelOf(el), el.qns.sym->GetPrincipleOffset()});
    }
    return rows;
}

std::vector<PolarizedEigenRow> PolarizedEigenRows(const EnergyLevels& combined, const EnergyLevels& up, const EnergyLevels& dn)
{
    // The HIGHEST occupied energy over BOTH channels: a doubly-empty level is dropped only when it sits ABOVE the
    // frontier, where it is one virtual among many.  BELOW the frontier it is a HOLE -- the whole point of looking.
    // (The old rule was silently inconsistent between cold and smeared runs: under Fermi smearing no occupation is
    // exactly 0.0, so nothing was ever dropped.)
    double eHomo=-1e300;
    for (const auto& elp : combined) if (elp.second.occ > 0.0) eHomo=std::max(eHomo, elp.second.e);
    std::set<Orbital_QNs> alreadyGotIt;
    std::vector<PolarizedEigenRow> rows;
    for (const auto& elp : combined)
    {
        const EnergyLevel& el=elp.second;
        Orbital_QNs upqns(el.qns.n, Spin::Up  , el.qns.sym);
        Orbital_QNs dnqns(el.qns.n, Spin::Down, el.qns.sym);
        if (alreadyGotIt.find(upqns)!=alreadyGotIt.end()) continue;
        alreadyGotIt.insert(upqns);
        alreadyGotIt.insert(dnqns);
        // A combined level need NOT exist in both spin channels (open shell).  Guard both lookups (find()==UB on a
        // miss in Release) and take the label + l from el.qns, whose sym is always valid.
        const EnergyLevel* u=up.FindOrNull(upqns);
        const EnergyLevel* d=dn.FindOrNull(dnqns);
        const double upOcc=u?u->occ:0.0, dnOcc=d?d->occ:0.0;
        if (upOcc==0.0 && dnOcc==0.0 && el.e>eHomo) continue;
        // ABSENT IS NOT EMPTY (2026-08-10, MnO run 29).  An absent channel used to fall back to the COMBINED level's
        // energy -- the OTHER channel's number -- and a fabricated occupancy of 0, so the row read as a level empty
        // AT THE SAME ENERGY with the splitting printing exactly 0.00000000; on MnO that manufactured a spin-up
        // "hole".  Now an absent level prints "--" and contributes no difference.
        const int upDeg=u?u->degen:(d?d->degen:1), dnDeg=d?d->degen:upDeg;
        PolarizedEigenRow r;
        r.occUp = u ? OccupationCell(upOcc, upDeg) : "--";
        r.eUp   = u ? Fixed8(u->e) : "--";
        r.label = LabelOf(el);
        r.occDn = d ? OccupationCell(dnOcc, dnDeg) : "--";
        r.eDn   = d ? Fixed8(d->e) : "--";
        r.dE    = (u && d) ? Fixed8(u->e-d->e) : "--";
        r.l     = el.qns.sym->GetPrincipleOffset();
        r.dnDim = (dnOcc==0.0);
        rows.push_back(std::move(r));
    }
    return rows;
}

} // namespace qchem::WaveFunction
