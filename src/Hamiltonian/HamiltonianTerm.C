// File: HamiltonianTermImplementation  General implementation of a HamiltonianTerm term in the Hamiltonian.


#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <ChargeDensity.H>
#include <Irrep_BS.H>
#include <Symmetry.H>
#include <iostream>
#include <cassert>

Static_HT_Imp::Static_HT_Imp()
    : itsExactCD(0)
{
    
};

Static_HT::SMat Static_HT_Imp::BuildHamiltonian(const TOrbital_IBS<double>* bs,const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    itsCache[qns]=CalculateHamiltonianMatrix(bs,s);
    itsBSs[qns]=bs;
 
    assert(itsCache.find(qns)!=itsCache.end());
    return itsCache[qns];
}

double Static_HT_Imp::CalculateEnergy(const Exact_CD* cd) const
{
    assert(itsExactCD);
    return itsExactCD->GetEnergy(this);
}

double Dynamic_HT_Imp::CalculateEnergy(const Exact_CD* cd) const
{
    if (true)
    {
        for (auto b:itsBSs)
            BuildHamiltonian(b.second,b.first.ms);
    }
    return itsExactCD->GetEnergy(this);
}

void Static_HT_Imp::UseChargeDensity(const Exact_CD* cd)
{
    itsExactCD =cd;
    assert(itsExactCD);
    
}
void Dynamic_HT_Imp::UseChargeDensity(const Exact_CD* cd)
{
    itsExactCD =cd;
    assert(itsExactCD);
    
}

const Static_HT_Imp::SMat& Static_HT_Imp::GetCachedMatrix(const TOrbital_IBS<double>* bs, const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    CacheMap::const_iterator i=itsCache.find(qns);
    assert(i!=itsCache.end());
    return i->second;
}

std::ostream&  Static_HT_Imp::Write(std::ostream& os) const
{
    return os;
}

std::istream&  Static_HT_Imp::Read (std::istream& is)
{
    return is;
}
