// File: HamiltonianTermImplementation  General implementation of a HamiltonianTerm term in the Hamiltonian.


#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <ChargeDensity.H>
#include <Irrep_BS.H>
#include <Symmetry.H>
#include <iostream>
#include <cassert>

Dynamic_HT_Imp::Dynamic_HT_Imp()
    : itsCD(0)
{
    
};

Static_HT::SMat Static_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    auto i=itsCache.find(qns);
    if (i==itsCache.end())
    {
        itsBSs[qns]=bs;    
        return itsCache[qns]=CalculateHamiltonianMatrix(bs,s);
    }
    else
        return i->second;
}

Dynamic_HT_Imp::SMat Dynamic_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    itsCache[qns]=CalculateHamiltonianMatrix(bs,s,itsCD);
    itsBSs[qns]=bs;
 
    assert(itsCache.find(qns)!=itsCache.end());
    return itsCache[qns];
}
Dynamic_HT_Imp::SMat Dynamic_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    itsCache[qns]=CalculateHamiltonianMatrix(bs,s,cd);
    itsBSs[qns]=bs;
 
    assert(itsCache.find(qns)!=itsCache.end());
    return itsCache[qns];
}

double Static_HT_Imp::CalculateEnergy(const DM_CD* cd) const
{
    return cd->GetEnergy(this);
}

double Dynamic_HT_Imp::CalculateEnergy(const DM_CD* cd) const
{
    assert(itsCD);
    if (true)
    {
        for (auto b:itsBSs)
            GetMatrix(b.second,b.first.ms);
    }
    return itsCD->GetEnergy(this);
}


void Dynamic_HT_Imp::UseChargeDensity(const DM_CD* cd)
{
    itsCD =cd;
    assert(itsCD);
    
}

bool Dynamic_HT_Imp::newCD(const DM_CD* cd) const
{
    if (cd==itsCD) 
        return false;
    else
    {
        itsCD=cd;
        return true;
    }

}


