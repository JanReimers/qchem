// File: Calculation/Imp/ValenceBasisGen.C  Implementation of the valence-basis generator (see ValenceBasisGen.C).
module;
#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cctype>

module qchem.ValenceBasisGen;

import qchem.PeriodicTable;                    // thePeriodicTable (symbol <-> Z)
import qchem.SCFIterator;                       // SCFParams (Verbose -> the per-iteration trace)
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH (default Zion)
import qchem.ChargeDensity;                    // Polarized_CD + Spin (the spin-resolved channel sampling)
import qchem.Types;                            // rvec3_t (radial sampling point)
import qchem.Math;                             // Pi (the 4 pi in the charge integral)

namespace qchem
{

std::vector<double> EvenTemperedWindow(int N, double emin, double emax)
{
    std::vector<double> es(N > 0 ? N : 0);
    const double beta = (N > 1) ? std::pow(emax/emin, 1.0/(N-1)) : 1.0;
    double e = emin;
    for (auto& a : es) { a = e; e *= beta; }
    return es;
}

// Append one uncontracted shell (one primitive per exponent) in Gaussian94 form.
static void WriteShell(std::ostringstream& os, int l, const std::vector<double>& es)
{
    const char L = "SPDFGH"[(l >= 0 && l < 6) ? l : 0];
    for (double e : es)
    {
        os << ' ' << L << "   1  1.00\n";
        os << "      " << std::fixed << std::setprecision(8) << std::setw(16) << e
           << "      " << std::setw(12) << 1.0 << '\n';
    }
}

GeneratedBasis GenerateValenceBasis(const ValenceBasisRecipe& r, bool showIterations)
{
    const int Z    = int(thePeriodicTable().GetZ(r.element));
    const int Zion = r.Zion > 0 ? r.Zion : Pseudopotential::GetGTH(r.element, r.functional).zion;
    const int nel  = r.electrons > 0 ? r.electrons : Zion;      // neutral pseudo-atom by default
    const int charge = Z - nel;                                 // AtomCalculation: electrons = Z - charge

    // --- validate: run the spherical pseudo-atom SCF in EXACTLY these per-l shells (LDA XC, GTH q=Zion).
    //     The atomic EC occupies only the l's the charge state fills; extra (polarization) l's ride along. ---
    AtomCalcOptions o;
    o.type            = AtomType::Gaussian;
    o.pseudopotential = true;
    o.valence         = Zion;
    o.exponentsByL    = r.shells;           // per-l independent lists: validate exactly what is emitted
    SCFParams p; p.Verbose = showIterations; // print the per-iteration convergence trace when asked
    p.MinVirial = 1e30;                       // the virial theorem does NOT hold under a pseudopotential (the PP
                                              //  replaces -Z/r), so |2+V/K| never -> 0: don't gate convergence on
                                              //  it here (this is the valence-gen backend only, NOT the system
                                              //  default -- the all-electron A_* path still uses the virial gate).
    AtomCalculation atom(Z, charge, o, p);

    // --- emit: the per-l shells as a Gaussian94 ELEMENT BLOCK (` EL 0 <shells> ****`) ---
    // The reader compares the file symbol against the UPPER-cased element (FindAtom), so upper-case it
    // (matching sipp.bsd's "SI"/"O").
    std::string sym = r.element;
    for (char& c : sym) c = char(std::toupper((unsigned char)c));
    std::ostringstream os;
    os << ' ' << sym << "   0\n";
    for (const auto& [l, es] : r.shells) WriteShell(os, l, es);
    os << " ****\n";

    return { atom.Energy(), atom.IsConverged(), os.str() };
}

double GenerateFloorEnergy(const ValenceBasisRecipe& r)
{
    const int Z    = int(thePeriodicTable().GetZ(r.element));
    const int Zion = r.Zion > 0 ? r.Zion : Pseudopotential::GetGTH(r.element, r.functional).zion;
    const int nel  = r.electrons > 0 ? r.electrons : Zion;
    const int charge = Z - nel;

    // Same element/PP/charge, but a large accuracy-POOL Gaussian basis (no exponentsByL) -- the near-complete
    // reference the recipe's window is measured against.
    AtomCalcOptions o;
    o.type            = AtomType::Gaussian;
    o.pseudopotential = true;
    o.valence         = Zion;
    o.accuracy        = BasisSetAccuracy::High;
    SCFParams p; p.MinVirial = 1e30;          // no virial gate under a PP (see GenerateValenceBasis)
    try   { AtomCalculation atom(Z, charge, o, p); return atom.Energy(); }
    catch (...) { return NAN; }               // pool conditioning can fail for some elements -> "floor unavailable"
}

GeneratedSeedDensity GenerateSeedDensity(const ValenceBasisRecipe& r, int Ngrid, double rmin, double rmax)
{
    const int Z    = int(thePeriodicTable().GetZ(r.element));
    const int Zion = r.Zion > 0 ? r.Zion : Pseudopotential::GetGTH(r.element, r.functional).zion;
    const int nel  = r.electrons > 0 ? r.electrons : Zion;      // valence electrons of THIS charge state
    const int charge = Z - nel;                                 // AtomCalculation: electrons = Z - charge

    // The SAME pseudo-atom SCF the basis generator runs (identical shells/PP/functional) -- one validated SCF.
    // spinResolved runs it POLARIZED: the atomic EC assigns the Hund unpaired electrons (e.g. Mn q7 -> S=5/2),
    // and the converged density is a Polarized_CD whose channels we sample alongside the total.
    AtomCalcOptions o;
    o.type            = AtomType::Gaussian;
    o.pseudopotential = true;
    o.valence         = Zion;
    o.exponentsByL    = r.shells;
    if (r.spinResolved) o.pol = Pol::Polarized;
    SCFParams p; p.MinVirial = 1e30;          // no virial gate under a PP (see GenerateValenceBasis)
    AtomCalculation atom(Z, charge, o, p);

    // The per-channel faces (spinResolved only): the polarized run's Density() IS-A Polarized_CD -- an
    // honest capability cross-cast (abstract face to abstract face, per the project cast rule).
    const ChargeDensity::Polarized_CD* pol = r.spinResolved
        ? dynamic_cast<const ChargeDensity::Polarized_CD*>(&atom.Density()) : nullptr;
    if (r.spinResolved && !pol)
        throw std::runtime_error("GenerateSeedDensity: polarized run did not yield a Polarized_CD for "+r.element);

    // Sample the spherical valence rho(r) (and the spin channels) on a LOG mesh (schema of
    // atomic_valence_densities.json), and integrate 4*pi*int r^2 rho dr (trapezoid) as the validation charge.
    const double beta = (Ngrid > 1) ? std::pow(rmax/rmin, 1.0/(Ngrid-1)) : 1.0;
    std::vector<double> rs(Ngrid), rho(Ngrid), up(Ngrid), dn(Ngrid);
    double rr = rmin;
    for (int i = 0; i < Ngrid; i++)
    {
        rs[i]  = rr;
        rho[i] = atom.Density()(rvec3_t(rr,0,0));
        if (pol)
        {
            up[i] = (*pol->GetChargeDensity(Spin::Up  ))(rvec3_t(rr,0,0));
            dn[i] = (*pol->GetChargeDensity(Spin::Down))(rvec3_t(rr,0,0));
        }
        rr *= beta;
    }
    auto radInt = [&](const std::vector<double>& f, int rpow)     // 4pi int r^(2+rpow) f dr, trapezoid
    {
        double sum = 0.0;
        for (int i = 0; i+1 < Ngrid; i++)
        {
            const double dr = rs[i+1]-rs[i];
            sum += 0.5*dr*(std::pow(rs[i],2+rpow)*f[i] + std::pow(rs[i+1],2+rpow)*f[i+1]);
        }
        return 4.0*Pi*sum;
    };
    const double q     = radInt(rho, 0);
    const double meanR = q > 0.0 ? radInt(rho, 1)/q : 0.0;                                // <r> (diffuseness)

    // UP-MAJORITY storage convention (doc/SCFSeedingPlan.md sec 10): rho_up is the MAJORITY channel; which
    // physical channel the SCF polarized is an arbitrary label, so swap if it landed the other way.
    double moment = 0.0;
    if (pol)
    {
        moment = radInt(up, 0) - radInt(dn, 0);                                           // ~ 2S
        if (moment < 0.0) { std::swap(up, dn); moment = -moment; }
    }

    // Emit the library entry (JSON object matching atomic_valence_densities.json).  rho stays the spin SUM
    // (spin-agnostic readers untouched); the pair rides as sibling arrays on the same grid.
    std::string sym = r.element;
    std::ostringstream os;
    os << std::setprecision(12);
    os << "{\"symbol\": \"" << sym << "\", \"Z\": " << Z << ", \"Nelec\": " << nel
       << ", \"charge\": " << q << ", \"functional\": \"" << r.functional << "\", \"kind\": \"valence\", "
       << "\"grid\": {\"kind\": \"log\", \"N\": " << Ngrid << ", \"rmin\": " << rmin << ", \"rmax\": " << rmax << "}, "
       << "\"rho\": [";
    for (int i = 0; i < Ngrid; i++) { if (i) os << ", "; os << rho[i]; }
    os << "]";
    if (pol)
    {
        os << ", \"moment\": " << moment;
        os << ", \"rho_up\": [";  for (int i = 0; i < Ngrid; i++) { if (i) os << ", "; os << up[i]; }  os << "]";
        os << ", \"rho_dn\": [";  for (int i = 0; i < Ngrid; i++) { if (i) os << ", "; os << dn[i]; }  os << "]";
    }
    os << "}";

    return { q, meanR, atom.IsConverged(), os.str(), moment };
}

std::string AssembleBasisFile(const std::string& name, const std::vector<std::string>& blocks)
{
    std::ostringstream os;
    os << name << "\n\n!\nBASIS=\"" << name << "\"\n";
    for (const std::string& b : blocks) os << b;
    return os.str();
}

} // namespace qchem
