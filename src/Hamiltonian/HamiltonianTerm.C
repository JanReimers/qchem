// File: HamiltonianTermImplementation  General implementation of a HamiltonianTerm term in the Hamiltonian.


#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <ChargeDensity.H>
#include <Symmetry.H>
#include <iostream>
#include <cassert>

HamiltonianTermImp::HamiltonianTermImp()
    : itsExactCD(0)
{
    MarkAllDirty();
};


void HamiltonianTermImp::MarkAllDirty()
{
    

}

HamiltonianTerm::SMat HamiltonianTermImp::BuildHamiltonian(const TOrbital_IBS<double>* bs,const Spin& s) const
{
    assert(bs);
    CacheIndex i(bs,s);
    itsCache[i]=CalculateHamiltonianMatrix(bs,s);
 
    assert(itsCache.find(i)!=itsCache.end());
    return itsCache[i];
}

double HamiltonianTermImp::CalculateEnergy() const
{
    if (DependsOnChargeDensity())
    {
        typedef CacheMap::const_iterator DITER;
        for (DITER i=itsCache.begin(); i!=itsCache.end(); i++)
        {
            
            const TOrbital_IBS<double>* bs=i->first.itsBasisSet;
            Spin s=i->first.itsSpin;
            assert(bs);
            BuildHamiltonian(bs,s);
            
        }
    }
    return itsExactCD->GetEnergy(this);
}

void HamiltonianTermImp::UseChargeDensity(const Exact_CD* cd)
{
//  TODO: SHould charge density have an unique ID?
//    assert(!itsExactCD || itsExactCD->GetID()!=theExactCD->GetID());
    itsExactCD =cd;
    assert(itsExactCD);
    if (DependsOnChargeDensity()) MarkAllDirty();
}

const HamiltonianTermImp::SMat& HamiltonianTermImp::GetCachedMatrix(const TOrbital_IBS<double>* bs, const Spin& s) const
{
    assert(bs);
    CacheIndex index(bs,s);
    CacheMap::const_iterator i=itsCache.find(index);
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
