// File: LDAVxc.C  Exact Coulomb potential



#include "HamiltonianImplementation/LDAVxc.H"
#include "BasisSet.H"
#include "oml/smatrix.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>

LDAVxc::LDAVxc()
    : HamiltonianTermImplementation( )
    , itsExchangeFunctional  (0)
{};

LDAVxc::LDAVxc(const rc_ptr<ExchangeFunctional>& lda)
    : HamiltonianTermImplementation(   )
    , itsExchangeFunctional  (lda)
{
    assert(&*itsExchangeFunctional);
};

void LDAVxc::UseChargeDensity(const ChargeDensity* exact)
{
    HamiltonianTermImplementation::UseChargeDensity(exact);
    itsExchangeFunctional->InsertChargeDensity(itsExactCD);
}

//########################################################################
//
//  Here Vxc is not fit to the exchange functional, so the Matrix and energy.
//  cannot be calculated analytically.
//
HamiltonianTerm::SMat LDAVxc::CalculateHamiltonianMatrix(const BasisSet* bs,const Spin&) const
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

void LDAVxc::InjectOverlaps(FittedFunction* ff, const BasisSet* theFitBasisSet) const
{
    const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(theFitBasisSet);
    assert(tbs);
    tbs->SetFitOverlap(ff,*itsExchangeFunctional);
}

void LDAVxc::InjectRepulsions(FittedFunction* ff, const BasisSet* theFitBasisSet) const
{
    const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(theFitBasisSet);
    assert(tbs);
    tbs->SetFitRepulsion(ff,*itsExchangeFunctional);
}

double LDAVxc::FitGetConstraint() const
{
    std::cerr << "There should be no constraint when fitting Vxc" << std::endl;
    exit(-1);
    return 0;
}


std::ostream& LDAVxc::Write(std::ostream& os) const
{
    os << *itsExchangeFunctional;

    return os;
}

std::istream& LDAVxc::Read (std::istream& is)
{
    itsExchangeFunctional.reset(ExchangeFunctional::Factory(is));
    is >> *itsExchangeFunctional;

    return is;
}

