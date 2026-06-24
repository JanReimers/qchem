// File: ElectronConfigurations/Imp/Crystal_EC.C  Electron configuration for a crystal (Bloch) calculation.
module;
#include <iostream>
module qchem.ElectronConfiguration.Crystal;

int Crystal_EC::GetN(const Irrep&) const {return itsNval;}

ElectronConfiguration::syms_t Crystal_EC::GetIrreps() const
{
    syms_t s;
    s.insert(itsIrrep.sym);
    return s;
}

void Crystal_EC::Display() const
{
    std::cout << "Crystal_EC: Nval=" << itsNval << std::endl;
}
