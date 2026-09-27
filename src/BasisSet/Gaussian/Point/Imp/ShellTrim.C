// File: BasisSet/Gaussian/Point/Imp/ShellTrim.C  The trim value and its reader decorator.
module;
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>
module qchem.BasisSet.Gaussian.Point.ShellTrim;
import qchem.PeriodicTable;   // the element symbol in the report
import qchem.Math;            // fabs
import qchem.Structure;       // Atom::itsZ

namespace qchem::BasisSet::Gaussian
{

bool ShellTrim::Removes(int Z, int l, const rvec_t& ex) const
{
    for (const auto& s : shells)
    {
        if (s.Z!=Z || s.l!=l || s.exponents.size()!=ex.size()) continue;
        bool same=true;
        for (size_t i=0;i<ex.size() && same;i++)
            same = fabs(s.exponents[i]-ex[i]) <= 1e-10*std::max(fabs(ex[i]),1.0);
        if (same) return true;
    }
    return false;
}

std::ostream& ShellTrim::Write(std::ostream& os) const
{
    if (shells.empty()) return os << "(none)";
    for (size_t i=0;i<shells.size();i++)
    {
        const auto& s=shells[i];
        os << (i?"; ":"") << thePeriodicTable().GetSymbol(s.Z) << " l=" << s.l << " {";
        for (size_t k=0;k<s.exponents.size();k++) os << (k?",":"") << s.exponents[k];
        os << "}";
    }
    return os;
}

// Read the next shell the TRIM keeps.  A multi-l shell (SP) survives with its untrimmed Ls; a shell whose
// every l is trimmed is skipped whole (the caller never sees it, so it is as if the file never held it).
GaussianRF* TrimmingReader::ReadNext(const Atom& atom)
{
    while (std::unique_ptr<GaussianRF> rf{itsInner.ReadNext(atom)})
    {
        const rvec_t ex=rf->GetExponents();
        itsLs.clear();
        for (int l : itsInner.GetLs()) if (!itsTrim.Removes(atom.itsZ, l, ex)) itsLs.push_back(l);
        if (!itsLs.empty()) return rf.release();   // the caller owns it, exactly as from the inner reader
    }
    itsLs.clear();
    return nullptr;
}

} // namespace
