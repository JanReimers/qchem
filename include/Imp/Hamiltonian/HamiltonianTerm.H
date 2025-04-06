// File: HamiltonianTermImplementation.H  General implementation of a HamiltonianTerm term in the Hamiltonian.
#ifndef _HamiltonianTermImp_H_
#define _HamiltonianTermImp_H_



#include <Hamiltonian.H>
#include <Spin.H>
#include <Orbital_QNs.H>

#include "oml/smatrix.h"
#include <map>

class HT_Common
{
protected:
    typedef std::map<Irrep_QNs,Static_HT::SMat> CacheMap;
    typedef std::map<Irrep_QNs,const Static_HT::ibs_t*> BSMap;
    mutable CacheMap   itsCache;       //Cache the H matricies for total energy calculations.
    mutable BSMap      itsBSs;
};


class Static_HT_Imp
    : public virtual Static_HT
    , protected HT_Common
{
public:
    virtual SMat GetMatrix(const ibs_t* bs,const Spin&) const;

protected:
    
    // Unconditional calculation, does no use cache.
    virtual SMat CalculateHamiltonianMatrix(const ibs_t*,const Spin&) const=0;
    // Calculation assumes that the cache is up to date. And E=Sum(Dab*Hab);
    virtual double CalculateEnergy(const DM_CD* cd) const;
   
  
};

class Dynamic_HT_Imp
: public virtual Dynamic_HT
, protected HT_Common
{
public:
    Dynamic_HT_Imp();
    virtual void UseChargeDensity(const DM_CD*);
    virtual SMat GetMatrix(const ibs_t* bs,const Spin&) const;
    virtual SMat GetMatrix(const ibs_t*,const Spin&,const DM_CD*) const; 

protected:
    // Unconditional calculation, does no use cache.
    virtual SMat CalculateHamiltonianMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const=0;
    virtual double CalculateEnergy(const DM_CD* cd) const;
    bool newCD(const DM_CD*) const;

    mutable const DM_CD*    itsCD;     //Density matrix charge density.
   
};



#endif
