// File: LDAVxc.C  Exact Coulomb potential



#include "Imp/Hamiltonian/LDAVxc.H"
#include "Imp/Fitting/FittedFunction.H"
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
    Dynamic_HT_Imp::UseChargeDensity(cd);
    itsExchangeFunctional->InsertChargeDensity(itsCD);
}

//########################################################################
//
//  Here Vxc is not fit to the exchange functional, so the Matrix and energy.
//  cannot be calculated analytically.
//
Static_HT::SMat LDAVxc::CalculateHamiltonianMatrix(const ibs_t* bs,const Spin&,const DM_CD* cd) const
{
    std::cerr << "LDAVxc::CalculateHamiltonianMatrix not implementated yet" << std::endl;
    exit(-1);
    return SMat();
}
void LDAVxc::GetEnergy(TotalEnergy&,const DM_CD* cd) const
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

