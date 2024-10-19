// File: FittedVee.h  Fitted Coulomb potential

#pragma implementation

#include "Imp/Hamiltonian/FittedVee.H"
#include "Imp/Hamiltonian/FittedVeeCD.H"
#include <TotalEnergy.H>
#include "oml/vector.h"
#include <cassert>

FittedVee::FittedVee()
    : PotentialImplementation( )
    , FittedFunctionImplementation<double>    ( )
{};

FittedVee::FittedVee(const rc_ptr<BasisSet>& bs)
    : PotentialImplementation(  )
    , FittedFunctionImplementation<double>(bs)
{};

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro(r2)/r12 .
//              /
//  Where ro is the charge density.
//
// The Hamiltonain matrix elements are calculated
//             /
//  Vee(i,j) = | dr Vfit(r) Oi(r) Oj(r) .
//             /
//
//           = Sum  { Ck <Oi|Vk|Oj> } .
//
//  This last part is carried out by the base class FitImplementation.
//
void FittedVee::BuildHamiltonian(const BasisSet* bs,const Spin&) const
{
    bs->AddOverlap(this);
}

void FittedVee::GetEnergy(TotalEnergy& te) const
{
    te.Eee = itsExactCD->GetOverlap(this);
}

void FittedVee::UseChargeDensity(const ChargeDensity& fitCD,const ChargeDensity& exactCD)
{
    PotentialImplementation::UseChargeDensity(fitCD,exactCD);
    FittablePotential* VeeTemp=new CDFittedVee;
    VeeTemp->UseChargeDensity(*itsFittedCD,*itsExactCD);
    DoFit(*VeeTemp);
    delete VeeTemp;
}

std::ostream& FittedVee::Write(std::ostream& os) const
{
    FittedFunctionImplementation<double>::Write(os);
    return os;
}

std::istream& FittedVee::Read (std::istream& is)
{
    FittedFunctionImplementation<double>::Read(is);
    return is;
}

Potential* FittedVee::Clone() const
{
    return new FittedVee(*this);
}

