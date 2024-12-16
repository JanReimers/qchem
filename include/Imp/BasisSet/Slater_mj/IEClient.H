// File: Slater_mj/IEClient.H
#ifndef _Slater_mj_IEClient_H_
#define _Slater_mj_IEClient_H_

#include "Imp/BasisSet/AtomIEClient.H"
    
namespace Slater_mj
{
 
struct IrrepIEClient 
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
{
    IrrepIEClient() {}; 
    IrrepIEClient(size_t N,int _kappa, double _mj);

    int kappa;
    double j,mj;
    
    void Init(const Vector  <double>& exponents);
    
private:
    virtual double Norm(double e, size_t l) const;
};

struct Dirac_IrrepIEClient
    : public virtual ::IrrepIEClient
{
    Dirac_IrrepIEClient() : itsLargeIEC(0), itsSmallIEC(0) {};
    virtual size_t size() const;

    void Init(const IrrepIEClient* large,const IrrepIEClient* small);
    
    const IrrepIEClient* itsLargeIEC;
    const IrrepIEClient* itsSmallIEC;
        
};


} //namespace Slater_mj

#endif // _Slater_mj_IEClient_H_
