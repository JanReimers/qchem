// File: Hamiltonian.H  Interface a Hamiltonianian operator.
#ifndef _Hamiltonian_H_
#define _Hamiltonian_H_

#include <types.H>
#include "Common/pmstream.h"
#include <ChargeDensity.H>
//
//  Abstract base for any HamiltonianTerm (HT) terms in the Hamiltonian.
//  We have two distinct types of HT:
//  Static_HT - Does not depend on the Charge Dnesity (CD), and therefore does change during iterations
//  Dynamic_HT - Vee, Vac which depend on the CD, and change with each iteration.
//
class Static_HT
    : public virtual PMStreamableObject
    , public virtual Static_CC
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;

    virtual const SMat& GetMatrix(const ibs_t*,const Spin&) const=0;
    virtual void        GetEnergy(EnergyBreakdown&,  const DM_CD*) const=0;
    virtual bool        IsPolarized() const {return false;}
};

class Dynamic_HT
    : public virtual PMStreamableObject
    , public virtual Dynamic_CC
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;    
    virtual const SMat& GetMatrix(const ibs_t*,const Spin&,const DM_CD*) const=0; 
    virtual void        GetEnergy(EnergyBreakdown&,  const DM_CD*) const=0;
    virtual bool        IsPolarized() const {return false;}
};

#include <Fitting/FittedFunctionClient.H>

class FittablePotential
    : public virtual Dynamic_HT
    , public virtual ScalarFFClient
{
public:
    virtual void UseChargeDensity(const DM_CD*)       =0;
};


class Hamiltonian
    : public virtual PMStreamableObject
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;

    virtual void            Add             (      Static_HT*)      =0;
    virtual void            Add             (      Dynamic_HT*)      =0;
    virtual SMat            GetMatrix(const ibs_t*,const Spin&,const DM_CD*)=0;
    virtual EnergyBreakdown GetTotalEnergy  (  const DM_CD*    ) const=0;
    virtual WaveFunction*   CreateWaveFunction(const BasisSet*,const ElectronConfiguration*,SCFAccelerator* acc) const=0;
    static  Hamiltonian*  Factory(std::istream&)        ;
};

#endif //_Hamiltonian_H_
