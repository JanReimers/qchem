// File: SphericalGaussianBS.H  Sphericalic gaussian basis set with orbital angular momentum.
#ifndef _SphericalGaussianBS_H_
#define _SphericalGaussianBS_H_



#include "IEClient.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "oml/vector.h"

namespace SphericalGaussian
{
   
class IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         TIrrepBasisSetCommon<double>
    , public         IrrepIEClient
{
public:
    IrrepBasisSet();
    IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,size_t size, double minexp, double maxexp, size_t L);
    IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Vector<double>& exponents, size_t L);

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    using IrrepBasisSetCommon::const_iterator;
    using IrrepBasisSetCommon::iterator;
    using IrrepBasisSetCommon::begin;
    using IrrepBasisSetCommon::end;


    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

} //namespace

#endif //_SphericalGaussianBS_H_
