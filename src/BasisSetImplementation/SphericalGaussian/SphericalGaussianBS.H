// File: SphericalGaussianBS.H  Sphericalic gaussian basis set with orbital angular momentum.
#ifndef _SphericalGaussianBS_H_
#define _SphericalGaussianBS_H_



#include "BasisSetImplementation/BasisSetImplementation.H"
#include "BasisSetImplementation/TBasisSetImplementation.H"

class SphericalGaussianBS
    : public virtual TBasisSet<double>
    , public         BasisSetImplementation
    , public         TBasisSetImplementation<double>
{
public:
    SphericalGaussianBS();
    SphericalGaussianBS(const LinearAlgebraParams&,IntegralDataBase<double>*,size_t size, double minexp, double maxexp, size_t L, Mesh* m=0);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual BasisSet* Clone(            ) const;
    virtual BasisSet* Clone(const RVec3&) const;
};

#endif //_SphericalGaussianBS_H_
