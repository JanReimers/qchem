// File: LDAVxc.C  Exact Coulomb potential



#include "Imp/Hamiltonian/LDAVxc.H"
#include "Imp/Fitting/FittedFunction.H"
#include "oml/smatrix.h"
#include <ChargeDensity.H>
#include <iostream>
#include <cassert>
#include <stdlib.h>

LDAVxc::LDAVxc()
    : HamiltonianTermImp( )
    , itsExchangeFunctional  (0)
{};

LDAVxc::LDAVxc(ex_t& lda)
    : HamiltonianTermImp(   )
    , itsExchangeFunctional  (lda)
{
    assert(&*itsExchangeFunctional);
};

void LDAVxc::UseChargeDensity(const Exact_CD* cd)
{
    HamiltonianTermImp::UseChargeDensity(cd);
    itsExchangeFunctional->InsertChargeDensity(itsExactCD);
}

//########################################################################
//
//  Here Vxc is not fit to the exchange functional, so the Matrix and energy.
//  cannot be calculated analytically.
//
HamiltonianTerm::SMat LDAVxc::CalculateHamiltonianMatrix(const TOrbital_IBS<double>* bs,const Spin&) const
{
    std::cerr << "LDAVxc::CalculateHamiltonianMatrix not implementated yet" << std::endl;
    exit(-1);
    return SMat();
}
void LDAVxc::GetEnergy(TotalEnergy&) const
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

