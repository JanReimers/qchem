// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


#include <AnalyticIE.H>
#include "oml/vector.h"



class SphericalGaussianIE1
    : public virtual AnalyticIE<double>
{

public:
    typedef AnalyticIE<double> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;
    typedef IE::iecv_t iecv_t;
    typedef Vector<int> IVec;
    
    struct SGparams
    {
        SGparams(const iecv_t&);
        size_t size() const {return Ls.size();}
        IVec Ls;  //Angular momentum.
        RVec es; //The orbital exponents.
        RVec ns; //Normalization constants
    };

    
        
    SphericalGaussianIE1() {};
    
    virtual size_t size() const {return -1;}

    virtual RVec MakeNormalization(iec_t* a) const;
    virtual RVec MakeCharge       (iec_t* a) const;
    
    virtual  SMat MakeOverlap  (iec_t* a         ) const;
    virtual  SMat MakeKinetic  (iec_t* a         ) const;
    virtual  SMat MakeNuclear  (iec_t* a,const Cluster&   ) const;
    virtual  SMat MakeRepulsion(iec_t* a         ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(iec_t*,iec_t*) const; //Used for charge and Vxc fitting.

    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT

    virtual void Make4C(ERI4& J, ERI4& K, const iecv_t&) const; //Used for Hartree-Fock

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

    virtual IE* Clone() const;

    typedef std::tuple<int,double,double> bf_tuple;
private:
    SMat MakeOverlap  (iec_t* ab,const bf_tuple& c) const;
    SMat MakeRepulsion(iec_t* ab,const bf_tuple& c) const;

    
    static double FourPi2;
};

#endif //_SphericalGaussianIE1_H_
