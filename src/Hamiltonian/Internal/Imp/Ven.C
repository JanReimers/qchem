// File: Ven.C  Electron-Nuclear potential.
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.ChargeDensity;

namespace qchem::Hamiltonian
{

Ven::Ven(const st_t& st)
    : rStatic_HT_Imp()
    , theStructure(st)
{
    assert(st->GetNumAtoms()>0);
};


rsmat_t Ven::MakeMatrix(const robs_t* bs,const Spin&) const
{
    // std::cout << "Ven=" << bs->Nuclear(&*theStructure) << std::endl;
    return bs->Nuclear(&*theStructure);
}

void Ven::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    const double e=cd->DM_Contract(this);
    te.Add("Een", e, EnergyRole::Potential, e);   // linear
}

std::ostream& Ven::Write(std::ostream& os) const
{
    os << "    Nuclear-electron potential Zi/|Ri-r| with " << theStructure->GetNumAtoms() << " atoms." << std::endl;
    return os;
}


} //namespace
