// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


#include "AnalyticIE.H"
#include "oml/vector.h"
#include <tuple>


class SphericalGaussianIE1
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
    typedef Vector<int> IVec;
    
    struct SGparams
    {
        SGparams(const iev_t&);
        size_t size() const {return Ls.size();}
        IVec Ls;  //Angular momentum.
        RVec es; //The orbital exponents.
        RVec ns; //Normalization constants
    };

    
        
    SphericalGaussianIE1() {};
    SphericalGaussianIE1(size_t L, const RVec& alphas);

    virtual size_t size() const {return es.size();}

   virtual RVec MakeNormalization() const;
    virtual RVec MakeCharge       () const;
    
    virtual  SMat MakeOverlap  (                 ) const;
    virtual  SMat MakeKinetic  (                 ) const;
    virtual  SMat MakeNuclear  (const Cluster&   ) const;
    virtual  SMat MakeRepulsion(                 ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(const AnalyticIE*) const; //Used for charge and Vxc fitting.

    virtual void MakeOverlap3C  (MList&,const AnalyticIE*) const; //Used for DFT
    virtual void MakeRepulsion3C(MList&,const AnalyticIE*) const; //Used for DFT
    virtual jk_t Make4C(const iev_t&) const; //Used for Hartree-Fock
    
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
