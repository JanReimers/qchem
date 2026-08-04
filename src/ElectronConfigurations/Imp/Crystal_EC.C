// File: ElectronConfigurations/Imp/Crystal_EC.C  Electron configuration for a crystal (Bloch) calculation.
module;
#include <iostream>
#include <vector>
module qchem.ElectronConfiguration.Crystal;

namespace qchem {

Crystal_EC::Crystal_EC(const Irrep& irr, int nval, bool globalFermi)
    : itsNval(nval), itsGlobalFermi(globalFermi)
{
    itsSyms.insert(irr.sym);
}

Crystal_EC::Crystal_EC(const std::vector<Irrep>& irreps, int nval, bool globalFermi)
    : itsNval(nval), itsGlobalFermi(globalFermi)
{
    for (const auto& irr : irreps) itsSyms.insert(irr.sym);
}

// ATOM SHELL CONVENTION (BlochQN doc, user design 2026-08-04): an IBZ wedge block carries its k-star
// multiplicity as its SPATIAL degeneracy (one stored representative per star), so in insulator mode the
// block fills star x Nval electrons -- every star copy's electrons flow through the one block, each band
// holding star x 2.  An unfolded mesh point (star = 1) is unchanged.  Metal (global-mu) mode instead asks
// for the WHOLE-MESH total (the composite solves one mu on the w x degeneracy-weighted capacities, which
// are convention-invariant), so it gets plain Nval regardless of which block asks.
int Crystal_EC::GetN(const Irrep& irr) const
{
    return itsGlobalFermi ? itsNval : itsNval*int(irr.sym->GetDegeneracy());
}

ElectronConfiguration::syms_t Crystal_EC::GetIrreps() const {return itsSyms;}

void Crystal_EC::Display() const
{
    std::cout << "Crystal_EC: Nval=" << itsNval
              << (itsGlobalFermi ? " total (global-μ metal), " : " per k-block, ")
              << itsSyms.size() << " k-block(s)" << std::endl;
}

} // namespace qchem