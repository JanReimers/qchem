// File: PolarizedGaussianIE.H  Integral Engine for Spherical gaussians.
#ifndef _PolarizedGaussianIE_H_
#define _PolarizedGaussianIE_H_

#include "RadialFunction.H"
#include "CDCache.H"
#include "Polarization.H"
#include <AnalyticIE.H>
#include "oml/vector.h"

namespace PolarizedGaussian
{

class IrrepIEClient;

class IntegralEngine
    : public virtual AnalyticIE<double>
{
public:
    typedef ScalarFunction <double> SF;
    typedef AnalyticIE<double>::Mat Mat;
    typedef AnalyticIE<double>::SMat SMat;
    typedef AnalyticIE<double>::RVec RVec;
    typedef const ::IrrepIEClient iec_t;

    IntegralEngine() {};
    virtual ~IntegralEngine() {};
    
    virtual RVec MakeNormalization(iec_t* a) const;
    virtual RVec MakeCharge       (iec_t* a) const;
    
    virtual  SMat MakeOverlap  (iec_t* a           ) const;
    virtual  SMat MakeKinetic  (iec_t* a         ) const;
    virtual  SMat MakeNuclear  (iec_t* a,const Cluster&   ) const;
    virtual  SMat MakeRepulsion(iec_t* a          ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(iec_t* a,iec_t* b) const; //Used for charge and Vxc fitting.

    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT
    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock
    
    virtual AnalyticIE<double>* Clone() const;
    virtual void Report(std::ostream&) const;
private:
    static const IrrepIEClient* dcast(iec_t*);
    SMat Integrate(qchem::IType2C , iec_t* a, const Cluster* cl=0) const;
    SMat Integrate(qchem::IType3C , iec_t* a, const RadialFunction* rc, const Polarization& pc) const;

    mutable CDCache cache; //Cache of all Gaussian pair charge distributions.
};
} //namespace PolarizedGaussian

#endif //_PolarizedGaussianIE_H_
