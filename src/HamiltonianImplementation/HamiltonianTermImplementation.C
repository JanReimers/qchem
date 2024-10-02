// File: HamiltonianTermImplementation  General implementation of a HamiltonianTerm term in the Hamiltonian.


#include "HamiltonianImplementation/HamiltonianTermImplementation.H"
#include "BasisSet.H"
#include "ChargeDensity.H"
#include "QuantumNumber.H"
#include <iostream>
#include <cassert>

HamiltonianTermImplementation::HamiltonianTermImplementation()
    : itsExactCD(0)
{
    MarkAllDirty();
};


void HamiltonianTermImplementation::MarkAllDirty()
{
    typedef DirtyMap::iterator DITER;
    for (DITER i=itsDirtyCache.begin(); i!=itsDirtyCache.end(); i++)
        i->second=true;

}

HamiltonianTerm::SMat HamiltonianTermImplementation::BuildHamiltonian(const BasisSet* bs,const Spin& s) const
{
    assert(bs);
    CacheIndex i(bs,s);
    bool notfound=itsDirtyCache.find(i)==itsDirtyCache.end();
    if (notfound || itsDirtyCache[i]==true)
    {
        itsCache[i]=CalculateHamiltonianMatrix(bs,s);
        itsDirtyCache[i]=false;
        StreamableObject::SetToPretty();
//        std::cout << "Recalculation H matrix for Qn=" << bs->GetQuantumNumber() << std::endl;
    }

    assert(itsCache.find(i)!=itsCache.end());
    return itsCache[i];
}

double HamiltonianTermImplementation::CalculateEnergy() const
{
    if (DependsOnChargeDensity())
    {
        typedef DirtyMap::const_iterator DITER;
        for (DITER i=itsDirtyCache.begin(); i!=itsDirtyCache.end(); i++)
        {
            if (i->second)
            {
                const BasisSet* bs=i->first.itsBasisSet;
                Spin s=i->first.itsSpin;
                assert(bs);
                BuildHamiltonian(bs,s);
            }
        }
    }
    return itsExactCD->GetEnergy(this);
}

void HamiltonianTermImplementation::UseChargeDensity(const ChargeDensity* theExactCD)
{
//  TODO: SHould charge density have an unique ID?
//    assert(!itsExactCD || itsExactCD->GetID()!=theExactCD->GetID());
    itsExactCD =theExactCD;
    assert(itsExactCD);
    if (DependsOnChargeDensity()) MarkAllDirty();
}

const HamiltonianTermImplementation::SMat& HamiltonianTermImplementation::GetCachedMatrix(const BasisSet* bs, const Spin& s) const
{
    assert(bs);
    CacheIndex index(bs,s);
    CacheMap::const_iterator i=itsCache.find(index);
    assert(i!=itsCache.end());
    return i->second;
}

std::ostream&  HamiltonianTermImplementation::Write(std::ostream& os) const
{
    return os;
}

std::istream&  HamiltonianTermImplementation::Read (std::istream& is)
{
    return is;
}
