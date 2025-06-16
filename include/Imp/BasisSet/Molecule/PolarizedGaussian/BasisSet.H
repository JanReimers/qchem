// File PolarizedGaussian/BasisSet.H
#ifndef _PolarizedGaussian_BasisSet_H_
#define _PolarizedGaussian_BasisSet_H_

#include "IEClient.H"
#include "../../BS_Common.H"
#include "DataBase/DB_HF.H"
#include "CDCache.H"

namespace PolarizedGaussian
{

class Reader;

class BasisSet 
    : public BS_Common
    , public DB_BS_2E<double>
{
public:
    BasisSet() {};
    BasisSet( Reader*, const Cluster*);
    virtual void Insert(bs_t* bs);


    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const;
private:
    mutable CDCache cache; //Cache of all Gaussian pair charge distributions.
};

} //namespace

#endif // _PolarizedGaussian

