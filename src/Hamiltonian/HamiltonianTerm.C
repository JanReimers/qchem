// File: HamiltonianTermImplementation  General implementation of a HamiltonianTerm term in the Hamiltonian.


#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <Irrep_BS.H>
#include <iostream>
#include <cassert>


Static_HT::SMat Static_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    auto i=itsCache.find(qns);
    if (i==itsCache.end())
    {
        itsBSs[qns]=bs;    
        return itsCache[qns]=CalcMatrix(bs,s);
    }
    else
        return i->second;
}

Dynamic_HT_Imp::SMat Dynamic_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    assert(bs);
    Irrep_QNs qns(s,&bs->GetQuantumNumber());
    auto i=itsCache.find(qns);
    if (itsCD!=cd || i==itsCache.end())
    {
        itsBSs[qns]=bs;  
        return itsCache[qns]=CalcMatrix(bs,s,cd);
    }
    else
        return i->second;
 
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


