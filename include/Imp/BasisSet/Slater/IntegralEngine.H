// File: SlaterIE.H  Integral Engine for Spherical gaussians.
#ifndef _SlaterIE_H_
#define _SlaterIE_H_


#include <AnalyticIE.H>
#include "Imp/BasisSet/AtomIE.H"
#include "oml/vector.h"


namespace Slater
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

            
    IntegralEngine();
    
    virtual size_t size() const {return -1;}

    //virtual  SMat MakeRepulsion(iec_t* a         ) const; //Used for charge and Vxc fitting.
    //virtual   Mat MakeRepulsion(iec_t*,iec_t*) const; //Used for charge and Vxc fitting.

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
    
    static double FourPi2;

};

} //namespace

#endif //_SlaterIE_H_
