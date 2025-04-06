// File: HamiltonianTermImplementation.H  General implementation of a HamiltonianTerm term in the Hamiltonian.
#ifndef _HamiltonianTermImp_H_
#define _HamiltonianTermImp_H_



#include <Hamiltonian.H>
#include <Spin.H>
#include <Orbital_QNs.H>

#include "oml/smatrix.h"
#include <map>

class Static_HT_Imp
    : public virtual Static_HT
{
public:
    virtual SMat GetMatrix(const ibs_t* bs,const Spin&) const;

protected:
    
    // Unconditional calculation, does no use cache.
    virtual SMat CalculateHamiltonianMatrix(const ibs_t*,const Spin&) const=0;
    // Calculation assumes that the cache is up to date. And E=Sum(Dab*Hab);
    virtual double CalculateEnergy(const DM_CD* cd) const;
   
    typedef std::map<Irrep_QNs,SMat> CacheMap;
    typedef std::map<Irrep_QNs,const ibs_t*> BSMap;
   

    mutable CacheMap   itsCache;       //Cache the H matricies for total energy calculations.
    mutable BSMap      itsBSs;
};

class Dynamic_HT_Imp
: public virtual Dynamic_HT
, public Static_HT_Imp
{
public:
    Dynamic_HT_Imp();
    virtual void UseChargeDensity(const DM_CD* exact);
    virtual SMat GetMatrix(const ibs_t* bs,const Spin&) const;
    
protected:
    virtual double CalculateEnergy(const DM_CD* cd) const;
    const DM_CD*    itsCD;     //Exact charge density.
   
};



#endif
