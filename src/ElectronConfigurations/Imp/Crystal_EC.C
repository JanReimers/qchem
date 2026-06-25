// File: ElectronConfigurations/Imp/Crystal_EC.C  Electron configuration for a crystal (Bloch) calculation.
module;
#include <iostream>
#include <vector>
module qchem.ElectronConfiguration.Crystal;

Crystal_EC::Crystal_EC(const Irrep& irr, int nval) : itsNval(nval)
{
    itsSyms.insert(irr.sym);
}

Crystal_EC::Crystal_EC(const std::vector<Irrep>& irreps, int nval) : itsNval(nval)
{
    for (const auto& irr : irreps) itsSyms.insert(irr.sym);
}

int Crystal_EC::GetN(const Irrep&) const {return itsNval;}

ElectronConfiguration::syms_t Crystal_EC::GetIrreps() const {return itsSyms;}

void Crystal_EC::Display() const
{
    std::cout << "Crystal_EC: Nval=" << itsNval << " per k-block, " << itsSyms.size() << " k-block(s)" << std::endl;
}
