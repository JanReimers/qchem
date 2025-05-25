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

const Static_HT::SMat& Static_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetSymmetry());
    auto i=itsCache.find(qns);
    if (i==itsCache.end())
        return itsCache[qns]=CalculateMatrix(bs,s);
    else
        return i->second;
}

const Dynamic_HT::SMat& Dynamic_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    assert(bs);
    
    Irrep_QNs qns(s,&bs->GetSymmetry());
    auto i=itsCache.find(qns);
    auto id=itsDirtyMap.find(qns);
    if (id==itsDirtyMap.end())
    {
        itsDirtyMap[qns]=true;
        id=itsDirtyMap.find(qns);
    }
    if (id->second || i==itsCache.end())
    {
        const SMat& H=itsCache[qns]=CalcMatrix(bs,s,cd); //This could reset itsDirtyMap[*]=false.
        itsDirtyMap[qns]=false; //Order it important for this line
        return H;
    }
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
        setDirty();
        return true;
    }

}

void Dynamic_HT_Imp::setDirty() const
{
    for (auto& d:itsDirtyMap) d.second=true;
}
const Dynamic_HT::SMat& Dynamic_HT_Imp_NoCache::GetMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    return itsMat=CalcMatrix(bs,s,cd);
}

