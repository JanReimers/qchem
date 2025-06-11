// File: HamiltonianTermImplementation.H  General implementation of a HamiltonianTerm term in the Hamiltonian.
#ifndef _HamiltonianTermImp_H_
#define _HamiltonianTermImp_H_



#include <Hamiltonian.H>
#include <Irrep_QNs.H>
#include "oml/smatrix.h"
#include <map>

class HT_Common
{
protected:
    typedef std::map<Irrep_QNs,Static_HT::SMat> CacheMap;
    typedef std::map<Irrep_QNs,const Static_HT::ibs_t*> BSMap;
    mutable CacheMap   itsCache;       //Cache the H matricies for total energy calculations.
};


class Static_HT_Imp
    : public virtual Static_HT
    , protected HT_Common
{
public:
    virtual const SMat& GetMatrix(const ibs_t* bs,const Spin&) const;

protected:
    // Unconditional calculation, does no use cache.
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const=0;
};

class Dynamic_HT_Imp
: public virtual Dynamic_HT
, protected HT_Common
{
public:
    Dynamic_HT_Imp();
    virtual const SMat& GetMatrix(const ibs_t*,const Spin&,const DM_CD*) const; 

protected:
    // Unconditional calculation, does not use cache.
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const=0;
    bool newCD(const DM_CD*) const;

    mutable const DM_CD* itsCD;      //Density matrix charge density.
};

// Used for polarized potentials (Vxc) which each polarization will handle its own cache.
class Dynamic_HT_Imp_NoCache
: public virtual Dynamic_HT
{
public:
    virtual const SMat& GetMatrix(const ibs_t*,const Spin&,const DM_CD*) const; 

protected:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const=0;
    mutable SMat itsMat;
};


#endif
