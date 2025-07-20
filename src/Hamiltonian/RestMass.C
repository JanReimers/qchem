// File: RestMass.C  Reast mass c^2 term for the Dirac hamiltonian.

#include <cassert>
#include <iostream>
#include <memory>
import qchem.DHF_IBS;
import qchem.ChargeDensity;
import qchem.Energy;
#include "RestMass.H"

import Common.Constants;

RestMass::RestMass()
    : Static_HT_Imp()
{};


Static_HT::SMat RestMass::CalculateMatrix(const ibs_t* bs,const Spin&) const
{
    static const double f=-2.0*c_light*c_light;
    // std::cout << "Rest mass/c^2=" << bs->GetRestMass() << std::endl;
    auto sbs=dynamic_cast<const Orbital_RKB_IBS<double>*>(bs);
    assert(sbs);
    return f*sbs->RestMass();
}

void RestMass::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    te.RestMass=cd->DM_Contract(this);
}

std::ostream& RestMass::Write(std::ostream& os) const
{
    os << "   Rest mass (beta-alpha)*c^2" << std::endl;
    return os;
}


