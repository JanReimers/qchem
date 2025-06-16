// File: LDAVxc.C  Exact Coulomb potential



#include "Imp/Hamiltonian/LDAVxc.H"
#include "oml/smatrix.h"
#include <ChargeDensity.H>
#include <iostream>
#include <cassert>
#include <stdlib.h>

LDAVxc::LDAVxc()
    : itsExchangeFunctional  (0)
{};

LDAVxc::LDAVxc(ex_t& lda)
    : itsExchangeFunctional  (lda)
{
    assert(&*itsExchangeFunctional);
};

void LDAVxc::UseChargeDensity(const DM_CD* cd)
{
    itsExchangeFunctional->InsertChargeDensity(cd);
}

//########################################################################
//
//  Here Vxc is not fit to the exchange functional, so the Matrix and energy.
//  cannot be calculated analytically.
//
Static_HT::SMat LDAVxc::CalcMatrix(const ibs_t* bs,const Spin&,const DM_CD* cd) const
{
    std::cerr << "LDAVxc::CalcMatrix not implementated yet" << std::endl;
    exit(-1);
    return SMat();
}
void LDAVxc::GetEnergy(EnergyBreakdown&,const DM_CD* cd) const
{
    std::cerr << "LDAVxc::GetEnergy(const BasisSet&,TotalEnergy&) not implementated yet" << std::endl;
    exit(-1);
}


std::ostream& LDAVxc::Write(std::ostream& os) const
{
    os << *itsExchangeFunctional;

    return os;
}

std::istream& LDAVxc::Read (std::istream& is)
{
    itsExchangeFunctional.reset(ExFunctional::Factory(is));
    is >> *itsExchangeFunctional;

    return is;
}

