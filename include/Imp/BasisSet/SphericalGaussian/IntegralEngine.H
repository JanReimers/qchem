// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


#include <AnalyticIE.H>
#include "Imp/BasisSet/AtomIE.H"
#include "oml/vector.h"


namespace SphericalGaussian
{
    
class IrrepIEClient;
    
class IntegralEngine
    : public virtual AnalyticIE<double>
    , public AtomIE
{

public:
    typedef AnalyticIE<double> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;

            
    IntegralEngine() {};
    
    virtual size_t size() const {return -1;}

    virtual  SMat MakeRepulsion(iec_t* a         ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(iec_t*,iec_t*) const; //Used for charge and Vxc fitting.

    //virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    //virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT

    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock

    virtual void Report(std::ostream&) const;
    virtual IE* Clone() const;

private:
    virtual double Overlap(double ea, double eb,size_t l) const;
    virtual double Kinetic(double ea ,double eb,size_t l) const;
    virtual double Nuclear(double ea, double eb,size_t l) const;
    virtual double Charge (double ea,           size_t l) const;
    virtual double Repulsion(double eab, double ec,size_t la,size_t lc) const;
    
    static const IrrepIEClient* dcast(iec_t*);
    typedef std::tuple<int,int,int,double,double> bf_tuple;
//    SMat MakeOverlap  (iec_t* ab,const bf_tuple& c) const;
    
    static double FourPi2;
};

} //namespace

#endif //_SphericalGaussianIE1_H_
