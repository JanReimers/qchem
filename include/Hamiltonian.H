// File: Hamiltonian.H  Interface a Hamiltonianian operator.
#ifndef _Hamiltonian_H_
#define _Hamiltonian_H_

#include "types.H"
#include "Imp/Misc/pmstream.h"


//
//  Abstract base for any HamiltonianTerm (HT) terms in the Hamiltonian.
//  We have two distinct types of HT:
//  Static_HT - Does not depend on the Charge Dnesity (CD), and therefore does change during iterations
//  Dynamic_HT - Vee, Vac which depend on the CD, and change with each iteration.
//
class Static_HT
    : public virtual PMStreamableObject
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;

    virtual SMat GetMatrix(const ibs_t*,const Spin&) const=0;
    virtual void GetEnergy(TotalEnergy&,  const DM_CD*) const=0;

    static  Static_HT*    Factory(std::istream&);

private:
    Static_HT& operator=(const Static_HT&);
};

class Dynamic_HT
: public virtual PMStreamableObject
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;    
    virtual SMat GetMatrix(const ibs_t*,const Spin&,const DM_CD*) const=0; 
    virtual void GetEnergy(TotalEnergy&,  const DM_CD*) const=0;
private:
    Dynamic_HT& operator=(const Dynamic_HT&);
};

#include <FittedFunctionClient.H>

class FittablePotential
    : public virtual Dynamic_HT
    , public virtual ScalarFFClient
{
public:
    virtual void UseChargeDensity(const DM_CD*)       =0;
    static  FittablePotential* Factory(std::istream&);
};


class Hamiltonian
    : public virtual PMStreamableObject
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;

    virtual void          Add             (      Static_HT*)      =0;
    virtual void          Add             (      Dynamic_HT*)      =0;
    virtual SMat          GetMatrix(const ibs_t*,const Spin&,const DM_CD*)=0;
    virtual TotalEnergy   GetTotalEnergy  (  const DM_CD*    ) const=0;
    
    static  Hamiltonian*  Factory(std::istream&)        ;
};

#endif //_Hamiltonian_H_
