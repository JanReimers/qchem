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
    typedef AnalyticIE<double> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;
    typedef IE::iecv_t iecv_t;
    typedef const ::IrrepIEClient iec_t;

    typedef std::vector<const Block*> blocks_t;
    
    IntegralEngine() {};
    IntegralEngine(const IrrepIEClient*);
    virtual ~IntegralEngine() {};
    
     struct PGparams
    {
        PGparams(const iecv_t&);
        size_t size() const {return ns.size();}
        std::vector<const RadialFunction*> radials; // Flattened radials
        std::vector<Polarization>          pols;    // Flattened polarizations
        Vector<double>                     ns;      //Norm constants
    };

    virtual RVec MakeNormalization(iec_t* a) const;
    virtual RVec MakeCharge       (iec_t* a) const;
    
    virtual  SMat MakeOverlap  (iec_t* a           ) const;
    virtual  SMat MakeKinetic  (iec_t* a         ) const;
    virtual  SMat MakeNuclear  (iec_t* a,const Cluster&   ) const;
    virtual  SMat MakeRepulsion(iec_t* a          ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(iec_t* a,iec_t* b) const; //Used for charge and Vxc fitting.

    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT
    virtual void Make4C(ERI4& J, ERI4& K, const iecv_t&) const; //Used for Hartree-Fock
    
    virtual IE* Clone() const;

private:
    static const IrrepIEClient* dcast(iec_t*);
    void Normalize(SMat&) const;
    SMat Integrate(RadialFunction::Types2C , iec_t* a, const Cluster* cl=0) const;
    SMat Integrate(RadialFunction::Types3C , iec_t* a, const RadialFunction* rc, const Polarization& pc) const;

    RVec                               ns;      //Norm constants
    SMat                               ons;     //Outer product of ns.
    mutable CDCache                    cache; //Cache of all Gaussian pair charge distributions.
};
} //namespace PolarizedGaussian

#endif //_PolarizedGaussianIE_H_
