// File: HamiltonianTermImplementation  General implementation of a HamiltonianTerm term in the Hamiltonian.


#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <ChargeDensity.H>
#include <Irrep_BS.H>
#include <Symmetry.H>
#include <iostream>
#include <cassert>

HamiltonianTermImp::HamiltonianTermImp()
    : itsExactCD(0)
{
    
};

HamiltonianTerm::SMat HamiltonianTermImp::BuildHamiltonian(const TOrbital_IBS<double>* bs,const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    itsCache[qns]=CalculateHamiltonianMatrix(bs,s);
    itsBSs[qns]=bs;
 
    assert(itsCache.find(qns)!=itsCache.end());
    return itsCache[qns];
}

double HamiltonianTermImp::CalculateEnergy() const
{
    if (DependsOnChargeDensity())
    {
        for (auto b:itsBSs)
            BuildHamiltonian(b.second,b.first.ms);
    }
    return itsExactCD->GetEnergy(this);
}

void HamiltonianTermImp::UseChargeDensity(const Exact_CD* cd)
{
//  TODO: SHould charge density have an unique ID?
//    assert(!itsExactCD || itsExactCD->GetID()!=theExactCD->GetID());
    itsExactCD =cd;
    assert(itsExactCD);
    
}

const HamiltonianTermImp::SMat& HamiltonianTermImp::GetCachedMatrix(const TOrbital_IBS<double>* bs, const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    CacheMap::const_iterator i=itsCache.find(qns);
    assert(i!=itsCache.end());
    return i->second;
}

std::ostream&  HamiltonianTermImp::Write(std::ostream& os) const
{
    return os;
}

std::istream&  HamiltonianTermImp::Read (std::istream& is)
{
    return is;
}
