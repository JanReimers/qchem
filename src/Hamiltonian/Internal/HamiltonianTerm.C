// File: HamiltonianTerm.C  General implementation of a HamiltonianTerm term in the Hamiltonian.
module;
#include <map>
export module qchem.Hamiltonian.Internal.Term;
export import qchem.Hamiltonian;
export import qchem.FittedFunctionClient;

import oml;
import qchem.Symmetry.Irrep;

export class HT_Common
{
protected:
    typedef std::map<Irrep_QNs,SMatrix<double>> CacheMap;
    // typedef std::map<Irrep_QNs,const Static_HT::ibs_t*> BSMap;
    mutable CacheMap   itsCache;       //Cache the H matrices for total energy calculations.
};


export class Static_HT_Imp
    : public virtual Static_HT
    , protected HT_Common
{
public:
    virtual const SMatrix<double>& GetMatrix(const ibs_t* bs,const Spin&) const;

protected:
    // Unconditional calculation, does no use cache.
    virtual SMatrix<double> CalculateMatrix(const ibs_t*,const Spin&) const=0;
};

export class Dynamic_HT_Imp
: public virtual Dynamic_HT
, protected HT_Common
{
public:
    Dynamic_HT_Imp();
    virtual const  SMatrix<double>& GetMatrix(const ibs_t*,const Spin&,const DM_CD*) const; 

protected:
    // Unconditional calculation, does not use cache.
    virtual SMatrix<double> CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const=0;
    bool newCD(const DM_CD*) const;

    mutable const DM_CD* itsCD;      //Density matrix charge density.
};

// Used for polarized potentials (Vxc) which each polarization will handle its own cache.
export class Dynamic_HT_Imp_NoCache
: public virtual Dynamic_HT
{
public:
    virtual const  SMatrix<double>& GetMatrix(const ibs_t*,const Spin&,const DM_CD*) const; 

protected:
    virtual SMatrix<double> CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const=0;
    mutable SMatrix<double> itsMat;
};

export class FittablePotential
    : public virtual Dynamic_HT
    , public virtual ScalarFFClient
{
public:
    virtual void UseChargeDensity(const DM_CD*)       =0;
};
