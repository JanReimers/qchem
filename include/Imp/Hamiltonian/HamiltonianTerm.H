// File: HamiltonianTermImplementation.H  General implementation of a HamiltonianTerm term in the Hamiltonian.
#ifndef _HamiltonianTermImp_H_
#define _HamiltonianTermImp_H_



#include <Hamiltonian.H>
#include <Spin.H>
#include <Orbital_QNs.H>

#include "oml/smatrix.h"
#include <map>

class HamiltonianTermImp
    : public virtual HamiltonianTerm
{
public:
    HamiltonianTermImp();

    virtual SMat BuildHamiltonian(const TOrbital_IBS<double>* bs,const Spin&) const;

    virtual void UseChargeDensity(const Exact_CD*);

    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

    const SMat& GetCachedMatrix(const TOrbital_IBS<double>* bs, const Spin& s) const;
protected:
    
    virtual bool DependsOnChargeDensity() const=0;
    // Unconditional calculation, does no use cache.
    virtual SMat CalculateHamiltonianMatrix(const TOrbital_IBS<double>*,const Spin&) const=0;
    // Calculation assumes that the cache is up to date. And E=Sum(Dab*Hab);
    double CalculateEnergy() const;
   
    typedef std::map<Irrep_QNs,SMat> CacheMap;
    typedef std::map<Irrep_QNs,const TOrbital_IBS<double>*> BSMap;
   

    const Exact_CD*    itsExactCD;     //Exact charge density.
    mutable CacheMap        itsCache;       //Cache the H matricies for total energy calculations.
    mutable BSMap     itsBSs;
};



#endif
