// File: PolarizedGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _PolarizedGaussianIE1_H_
#define _PolarizedGaussianIE1_H_

#include "BasisSetImplementation/PolarizedGaussian/RadialFunction.H"
#include "AnalyticIE.H"
#include "Misc/Polarization.H"
#include "oml/vector.h"


class BasisFunctionBlock;

class PolarizedGaussianIE1
    : public virtual AnalyticIE<double>
{
public:
    typedef AnalyticIE<double> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;
    typedef IE::iev_t iev_t;
    typedef IE::jk_t jk_t;
    typedef const IrrepIEClient iec_t;

    typedef std::vector<const BasisFunctionBlock*> blocks_t;
    
    PolarizedGaussianIE1() {};
    PolarizedGaussianIE1(const blocks_t&);
    virtual ~PolarizedGaussianIE1() {};

    virtual size_t size() const {return radials.size();}

    virtual RVec MakeNormalization() const;
    virtual RVec MakeCharge       (iec_t* a) const;
    
    virtual  SMat MakeOverlap  (iec_t* a           ) const;
    virtual  SMat MakeKinetic  (iec_t* a         ) const;
    virtual  SMat MakeNuclear  (iec_t* a,const Cluster&   ) const;
    virtual  SMat MakeRepulsion(iec_t* a          ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(const IE*a,const IE*b) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(iec_t* a,iec_t* b) const; //Used for charge and Vxc fitting.

    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT
    virtual jk_t Make4C(const iev_t&) const; //Used for Hartree-Fock
    
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

    virtual IE* Clone() const;

private:
    void Normalize(SMat&) const;
    SMat Integrate(RadialFunction::Types2C , iec_t* a, const Cluster* cl=0) const;
    SMat Integrate(RadialFunction::Types3C , iec_t* a, const RadialFunction* rc, const Polarization& pc) const;

    SMat Integrate(RadialFunction::Types2C , const Cluster* cl=0) const;
    SMat Integrate(RadialFunction::Types3C , const RadialFunction* rc, const Polarization& pc) const;

    blocks_t                           blocks;  // by block one radial function for many polarizations
    std::vector<const RadialFunction*> radials; // Flattened radials
    std::vector<Polarization>          pols;    // Flattened polarizations
    RVec                               ns;      //Norm constants
    SMat                               ons;     //Outer product of ns.
    mutable CDcache                    cache; //Cache of all Gaussian pair charge distributions.
};

#endif //_PolarizedGaussianIE1_H_
