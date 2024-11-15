// File: SphericalGaussian_mIE.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussian_mIE_H_
#define _SphericalGaussian_mIE_H_

#include "Imp/BasisSet/AtomIE.H"
#include <AnalyticIE.H>
#include "oml/vector.h"


namespace SphericalGaussian_m
{
    
class IrrepIEClient;
class IEClient;
    
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

    virtual RVec MakeNormalization(iec_t* a) const;
    virtual RVec MakeCharge       (iec_t* a) const;
    
    //virtual  SMat MakeOverlap  (iec_t* a         ) const;
    //virtual  SMat MakeKinetic  (iec_t* a         ) const;
    virtual  SMat MakeNuclear  (iec_t* a,const Cluster&   ) const;
    virtual  SMat MakeRepulsion(iec_t* a         ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(iec_t*,iec_t*) const; //Used for charge and Vxc fitting.

    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT

    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock
    virtual void MakeDirect  (erij_t& Jab, const ::IEClient*) const; //Used for Hartree-Fock
    virtual void MakeExchange(erik_t& Kab, const ::IEClient*) const; //Used for Hartree-Fock

    virtual void Report(std::ostream&) const;
    virtual IE* Clone() const;

private:
    virtual double Overlap(double eab,size_t lab) const;
    virtual double Kinetic(double ea, double eb,size_t l) const;

    
    static const IrrepIEClient* dcast(iec_t*);
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    SMat MakeOverlap  (iec_t* ab,const bf_tuple& c) const;
    SMat MakeRepulsion(iec_t* ab,const bf_tuple& c) const;
    ERIJ MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c,const IEClient*) const;
    ERIK MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b,const IEClient*) const;
   
    static double FourPi2;
};

} //namespace

#endif //_SphericalGaussian_mIE_H_
