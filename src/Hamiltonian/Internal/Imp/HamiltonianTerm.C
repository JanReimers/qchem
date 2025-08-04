// File: HamiltonianTerm  General implementation of a HamiltonianTerm term in the Hamiltonian.
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
module qchem.Hamiltonian.Internal.Term;
import qchem.ChargeDensity;
import qchem.IrrepBasisSet;
import qchem.Symmetry;

Dynamic_HT_Imp::Dynamic_HT_Imp() : itsCD(0) {};

const SMatrix<double>& Static_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s) const
{
    assert(bs);
    Irrep_QNs qns(s,bs->GetSymmetry());
    auto i=itsCache.find(qns);
    if (i==itsCache.end())
        return itsCache[qns]=CalculateMatrix(bs,s);
    else
        return i->second;
}

const SMatrix<double>& Dynamic_HT_Imp::GetMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
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

const SMatrix<double>& Dynamic_HT_Imp_NoCache::GetMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    return itsMat=CalcMatrix(bs,s,cd);
}

