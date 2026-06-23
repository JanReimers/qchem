// File: HamiltonianTerm.C  General implementation of a HamiltonianTerm term in the Hamiltonian.
module;
#include <map>
export module qchem.Hamiltonian.Internal.Term;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;

import qchem.Symmetry.Irrep;

export namespace qchem::Hamiltonian
{

class HT_Common
{
protected:
    typedef std::map<Irrep,rsmat_t> CacheMap;
    // typedef std::map<Irrep,const Static_HT::obs_t*> BSMap;
    mutable CacheMap   itsCache;       //Cache the H matrices for total energy calculations.
};


class Static_HT_Imp
    : public virtual Static_HT
    , protected HT_Common
{
public:
    virtual const rsmat_t& GetMatrix(const obs_t* bs,const Spin&) const;

protected:
    // Unconditional calculation, does no use cache.
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const=0;
};

class Dynamic_HT_Imp
: public virtual Dynamic_HT
, protected HT_Common
{
public:
    Dynamic_HT_Imp();
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const DM_CD*) const; 

protected:
    // Unconditional calculation, does not use cache.
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const DM_CD* cd) const=0;
    bool newCD(const DM_CD*) const;

    mutable const DM_CD* itsCD;      //Density matrix charge density.
};

// Used for polarized potentials (Vxc) which each polarization will handle its own cache.
class Dynamic_HT_Imp_NoCache
: public virtual Dynamic_HT
{
public:
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const DM_CD*) const; 

protected:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const DM_CD* cd) const=0;
    mutable rsmat_t itsMat;
};

// A density-dependent Hamiltonian potential term.  (No longer a Fitting::ScalarFFClient -- the
// "I can answer the fitter's questions" role belongs on the concrete LDAVxc, not this term interface.)
class FittablePotential
    : public virtual Dynamic_HT
{
public:
    virtual void UseChargeDensity(const DM_CD*)       =0;
};

} //namespace
