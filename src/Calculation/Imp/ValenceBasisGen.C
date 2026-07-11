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
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH (default Zion)

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

GeneratedBasis GenerateValenceBasis(const ValenceBasisRecipe& r)
{
    const int Z    = int(thePeriodicTable().GetZ(r.element));
    const int Zion = r.Zion > 0 ? r.Zion : Pseudopotential::GetGTH(r.element, r.functional).zion;
    const int nel  = r.electrons > 0 ? r.electrons : Zion;      // neutral pseudo-atom by default
    const int charge = Z - nel;                                 // AtomCalculation: electrons = Z - charge

    // --- validate: run the spherical pseudo-atom SCF in the validation window (LDA XC, GTH q=Zion) ---
    AtomCalcOptions o;
    o.type            = AtomType::Gaussian;
    o.pseudopotential = true;
    o.valence         = Zion;
    o.exponents       = r.validateWindow;   // "bring your own exponents" (ltrim=0 -> every occupied l gets it)
    AtomCalculation atom(Z, charge, o);

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

std::string AssembleBasisFile(const std::string& name, const std::vector<std::string>& blocks)
{
    std::ostringstream os;
    os << name << "\n\n!\nBASIS=\"" << name << "\"\n";
    for (const std::string& b : blocks) os << b;
    return os.str();
}

} // namespace qchem
