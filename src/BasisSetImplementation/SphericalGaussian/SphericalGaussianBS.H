// File: SphericalGaussianBS.H  Sphericalic gaussian basis set with orbital angular momentum.
#ifndef _SphericalGaussianBS_H_
#define _SphericalGaussianBS_H_



#include "BasisSetImplementation/BasisSetImplementation.H"
#include "BasisSetImplementation/TBasisSetImplementation.H"
#include "oml/vector.h"

struct SphericalGaussianIEClient
    : public virtual IrrepIEClient
{
    SphericalGaussianIEClient() {};
    SphericalGaussianIEClient(size_t N, double minexp,double maxexp,size_t L);
    Vector<size_t> Ls;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
};
    

class SphericalGaussianBS
    : public virtual TIrrepBasisSet<double>
    , public         BasisSetImplementation
    , public         TBasisSetImplementation<double>
    , public SphericalGaussianIEClient
{
public:
    SphericalGaussianBS();
    SphericalGaussianBS(const LinearAlgebraParams&,IntegralDataBase<double>*,size_t size, double minexp, double maxexp, size_t L, Mesh* m=0);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual IrrepBasisSet* Clone(            ) const;
    virtual IrrepBasisSet* Clone(const RVec3&) const;
};

#endif //_SphericalGaussianBS_H_
