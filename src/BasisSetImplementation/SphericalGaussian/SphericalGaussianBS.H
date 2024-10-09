// File: SphericalGaussianBS.H  Sphericalic gaussian basis set with orbital angular momentum.
#ifndef _SphericalGaussianBS_H_
#define _SphericalGaussianBS_H_



#include "BasisSetImplementation/SphericalGaussian/IEClient.H"
#include "BasisSetImplementation/BasisSetImplementation.H"
#include "BasisSetImplementation/TBasisSetImplementation.H"
#include "oml/vector.h"


class SphericalGaussianBS
    : public virtual TIrrepBasisSet<double>
    , public         BasisSetImplementation
    , public         TBasisSetImplementation<double>
    , public SphericalGaussianIEClient
{
public:
    SphericalGaussianBS();
    SphericalGaussianBS(const LinearAlgebraParams&,IntegralDataBase<double>*,size_t size, double minexp, double maxexp, size_t L);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual IrrepBasisSet* Clone(            ) const;
    virtual IrrepBasisSet* Clone(const RVec3&) const;
};

#endif //_SphericalGaussianBS_H_
