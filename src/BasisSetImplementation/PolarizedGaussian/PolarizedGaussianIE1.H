// File: PolarizedGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _PolarizedGaussianIE1_H_
#define _PolarizedGaussianIE1_H_

#include "BasisSetImplementation/PolarizedGaussian/RadialFunction.H"
#include "IntegralEngine1.H"
#include "Misc/Polarization.H"
#include "oml/vector.h"


class BasisFunctionBlock;

class PolarizedGaussianIE1
    : public virtual IntegralEngine1<double>
{
public:
    typedef IntegralEngine1<double> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;
    typedef IE::iev_t iev_t;
    typedef IE::jk_t jk_t;
    typedef std::vector<const BasisFunctionBlock*> blocks_t;
    
    PolarizedGaussianIE1() {};
    PolarizedGaussianIE1(const blocks_t&, const RVec& ns);
    virtual ~PolarizedGaussianIE1() {};

    virtual size_t size() const {return ns.size();}

    virtual SMat MakeOverlap    (         ) const;
    virtual Mat  MakeOverlap    (const IE*) const;
    virtual RVec MakeOverlap    (const SF&) const;
    virtual void MakeOverlap3C  (MList&,const IE*) const;
//
    virtual SMat MakeRepulsion  (          ) const;
    virtual Mat  MakeRepulsion  (const IE* ) const;
    virtual Vector<double>  MakeRepulsion  (const SF&  ) const;
    virtual void MakeRepulsion3C(MList&,const IE*) const;
    
    virtual void MakeRepulsion4C(ERIList& Coulomb, ERIList& exchange, const iev_t&) const;
    virtual jk_t Make4C(const iev_t&) const;

    virtual SMat MakeKinetic  (              ) const;
    virtual SMat MakeNuclear  (const Cluster&) const;
    virtual RVec MakeNormalization() const;
    virtual RVec MakeCharge       () const;
    
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

    virtual IE* Clone() const;

private:
    void Normalize(SMat&) const;
    SMat MakeOverlap3C  (const RadialFunction* rc, const Polarization& pc) const;
    SMat MakeRepulsion3C(const RadialFunction* rc, const Polarization& pc) const;


    blocks_t                           blocks;  // by block one radial function for many polarizations
    std::vector<const RadialFunction*> radials; // Flattened radials
    std::vector<Polarization>          pols;    // Flattened polarizations
    RVec                               ns;      //Norm constants
    SMat                               ons;     //Outer product of ns.
    mutable CDcache                    cache; //Cache of all Gaussian pair charge distributions.
};

#endif //_PolarizedGaussianIE1_H_
