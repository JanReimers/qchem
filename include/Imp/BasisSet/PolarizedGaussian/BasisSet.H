// File PolarizedGaussian/BasisSet.H
#ifndef _PolarizedGaussian_BasisSet_H_
#define _PolarizedGaussian_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"
#include "Imp/DataBase/DB_HF.H"
#include "CDCache.H"

namespace PolarizedGaussian
{

class Reader;

class BasisSet 
    : public BasisSetImp
    , public DB_BS_2E<double>
{
public:
    BasisSet() {};
    BasisSet( Reader*, const Cluster*);
    //! Make a primatives dentered at R=(0.0,0) with polarizations up to LMax
    BasisSet(size_t N, double emin, double emax, size_t LMax, const Cluster*);
    virtual void Insert(bs_t* bs);


    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const;
private:
    mutable CDCache cache; //Cache of all Gaussian pair charge distributions.
};

} //namespace

#endif // _PolarizedGaussian

