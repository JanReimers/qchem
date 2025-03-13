// File: PolarizedGaussianIE.H  Integral Engine for Spherical gaussians.
#ifndef _PolarizedGaussianIE_H_
#define _PolarizedGaussianIE_H_

#include "RadialFunction.H"
#include "CDCache.H"
#include "Polarization.H"
#include "Imp/DataBase/HeapDB.H"

#include <AnalyticIE.H>
#include "oml/vector.h"

namespace PolarizedGaussian
{

class IrrepIEClient;

class IntegralEngine1
    : public virtual Integrals_1E<double>
    , public DB_1E<double>
    , public DB_2E<double>
    , public DB_DFT<double>
 
{
public:
protected:
    IntegralEngine1() {};
    IntegralEngine1(const DB_BS_2E<double>* db) : DB_2E<double>(db) {};
protected:
    virtual Vec  MakeCharge() const;


    typedef Integrals_Base<double> Base;
    //typedef DB_Fit::Vec Vec;
    typedef Base::SMat SMat;
    typedef Base::obs_t obs_t;
    typedef Base::fbs_t fbs_t; //Fit basis
    typedef Base::Mat Mat;
    typedef Integrals_DFT<double>::ERI3 ERI3;

    virtual SMat MakeOverlap() const {return MakeIntegrals(qchem::Overlap2C);}
    virtual SMat MakeRepulsion() const {return MakeIntegrals(qchem::Repulsion2C);}
    virtual SMat MakeKinetic() const {return MakeIntegrals(qchem::Kinetic);}
    virtual SMat MakeNuclear(const Cluster* cl) const {return MakeIntegrals(qchem::Nuclear,cl);}
    virtual Mat  MakeRepulsion(const fbs_t& b) const;
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const; //Used for DFT
    virtual ERI4 MakeDirect  (const obs_t& c) const;
    virtual ERI4 MakeExchange(const obs_t& b) const;

    SMat MakeIntegrals(qchem::IType2C,const Cluster*cl =0) const;
    SMat Integrate(qchem::IType3C , const RadialFunction* rc, const Polarization& pc) const;
    mutable CDCache cache; //Cache of all Gaussian pair charge distributions.
};

} //namespace PolarizedGaussian

#endif //_PolarizedGaussianIE_H_
