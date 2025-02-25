// File: HamiltonianTermImplementation  General implementation of a HamiltonianTerm term in the Hamiltonian.


#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <BasisSet.H>
#include <ChargeDensity.H>
#include <QuantumNumber.H>
#include <iostream>
#include <cassert>

HamiltonianTermImp::HamiltonianTermImp()
    : itsExactCD(0)
{
    MarkAllDirty();
};


void HamiltonianTermImp::MarkAllDirty()
{
    typedef DirtyMap::iterator DITER;
    for (DITER i=itsDirtyCache.begin(); i!=itsDirtyCache.end(); i++)
        i->second=true;

}

HamiltonianTerm::SMat HamiltonianTermImp::BuildHamiltonian(const Orbital_IBS<double>* bs,const Spin& s) const
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

double HamiltonianTermImp::CalculateEnergy() const
{
    if (DependsOnChargeDensity())
    {
        typedef DirtyMap::const_iterator DITER;
        for (DITER i=itsDirtyCache.begin(); i!=itsDirtyCache.end(); i++)
        {
            if (i->second)
            {
                const Orbital_IBS<double>* bs=i->first.itsBasisSet;
                Spin s=i->first.itsSpin;
                assert(bs);
                BuildHamiltonian(bs,s);
            }
        }
    }
    return itsExactCD->GetEnergy(this);
}

void HamiltonianTermImp::UseChargeDensity(const ChargeDensity* theExactCD)
{
//  TODO: SHould charge density have an unique ID?
//    assert(!itsExactCD || itsExactCD->GetID()!=theExactCD->GetID());
    itsExactCD =theExactCD;
    assert(itsExactCD);
    if (DependsOnChargeDensity()) MarkAllDirty();
}

const HamiltonianTermImp::SMat& HamiltonianTermImp::GetCachedMatrix(const Orbital_IBS<double>* bs, const Spin& s) const
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
