// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


#include "IntegralEngine1.H"
#include "oml/vector.h"
#include <tuple>

class SphericalGaussianIE1
    : public virtual IntegralEngine1<double>
{
public:
    typedef IntegralEngine1<double> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;
    typedef IE::iev_t iev_t;
    
    SphericalGaussianIE1(size_t L, const RVec& alphas);

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
//
    virtual SMat MakeKinetic  (              ) const;
    virtual SMat MakeNuclear  (const Cluster&) const;
    virtual RVec MakeNormalization() const;
    virtual RVec MakeCharge       () const;
    
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

    virtual IE* Clone() const;

    typedef std::tuple<int,double,double> bf_tuple;
    bf_tuple operator()(index_t i) const {return std::make_tuple(L,es(i),ns(i));}
private:
    SMat MakeOverlap  (const bf_tuple&) const;
    SMat MakeRepulsion(const bf_tuple&) const;

    int  L;  //Angular momentum.
    RVec es; //The orbital exponents.
    RVec ns; //Normalization constants
   
    static double FourPi2;
};

#endif //_SphericalGaussianIE1_H_
