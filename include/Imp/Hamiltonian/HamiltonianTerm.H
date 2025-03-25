// File: HamiltonianTermImplementation.H  General implementation of a HamiltonianTerm term in the Hamiltonian.
#ifndef _HamiltonianTermImp_H_
#define _HamiltonianTermImp_H_



#include <Hamiltonian.H>
#include <Spin.H>
#include "oml/smatrix.h"
#include <map>

class HamiltonianTermImp
    : public virtual HamiltonianTerm
{
public:
    HamiltonianTermImp();

    virtual SMat BuildHamiltonian(const TOrbital_IBS<double>* bs,const Spin&) const;

    virtual void UseChargeDensity(const Exact_CD* exact);

    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

    const SMat& GetCachedMatrix(const TOrbital_IBS<double>* bs, const Spin& s) const;
protected:
    class CacheIndex
    {
    public:
        CacheIndex(const TOrbital_IBS<double>* bs, const Spin& s)
        : itsBasisSet(bs)
        , itsSpin(s)
        {};


        friend bool operator<(const CacheIndex& a, const CacheIndex& b)
        {
            if (a.itsBasisSet == b.itsBasisSet )
            {
                return  a.itsSpin < b.itsSpin;
            }
            else
            {
                return a.itsBasisSet < b.itsBasisSet;
            }
        }

        const TOrbital_IBS<double>* itsBasisSet;
        Spin            itsSpin;
    };
    virtual bool DependsOnChargeDensity() const=0;
    // Unconditional calculation, does no use cache.
    virtual SMat CalculateHamiltonianMatrix(const TOrbital_IBS<double>*,const Spin&) const=0;
    // Calculation assumes that the cache is up to date. And E=Sum(Dab*Hab);
    double CalculateEnergy() const;
    void MarkAllDirty();

    typedef std::map<CacheIndex,SMat> CacheMap;
    typedef std::map<CacheIndex,bool> DirtyMap;

    const Exact_CD*    itsExactCD;     //Exact charge density.
    mutable CacheMap        itsCache;       //Cache the H matricies for total energy calculations.
    mutable DirtyMap        itsDirtyCache;  //Is cache is older than charge density.
};



#endif
