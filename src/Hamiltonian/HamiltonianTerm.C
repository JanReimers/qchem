// File: HamiltonianTermImplementation  General implementation of a HamiltonianTerm term in the Hamiltonian.


#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <ChargeDensity.H>
#include <BasisSet/Irrep_BS.H>
#include <Symmetry.H>
#include <iostream>
#include <cassert>

Dynamic_HT_Imp::Dynamic_HT_Imp()
    : itsCD(0)
{
    
};

const Static_HT::SMat& Static_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,bs->GetSymmetry());
    auto i=itsCache.find(qns);
    if (i==itsCache.end())
        return itsCache[qns]=CalculateMatrix(bs,s);
    else
        return i->second;
}

const Dynamic_HT::SMat& Dynamic_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    assert(bs);
    Irrep_QNs qns(s,bs->GetSymmetry());
    if (auto i=itsCache.find(qns);i==itsCache.end())
        return itsCache[qns]=CalcMatrix(bs,s,cd); //This could cler the cache if cd is new.
    else
        return i->second; //Cache version 
}

bool Dynamic_HT_Imp::newCD(const DM_CD* cd) const
{
    assert(cd);
    if (cd==itsCD) 
        return false;
    else
    {
        itsCD=cd;
        itsCache.clear();
        return true;
    }

}

const Dynamic_HT::SMat& Dynamic_HT_Imp_NoCache::GetMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    return itsMat=CalcMatrix(bs,s,cd);
}

