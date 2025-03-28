// File: Atom/ml/Gaussian_BF.H r^l exp(-a*r^2) type Gaussian basis function.
#ifndef _Atom_ml_Gaussian_BF_H_
#define _Atom_ml_Gaussian_BF_H_

#include "../l/Gaussian_BF.H"

namespace Atom_ml
{
namespace Gaussian
{

    class BasisFunction : public Atoml::Gaussian::BasisFunction
    {
    public:
        BasisFunction(double theExponent,int n, int l, int m, double norm);
        
        virtual bool   operator==(const ::BasisFunction&) const;
    
        virtual BasisFunction* Clone(        ) const;
    private:
        typedef Atoml::Gaussian::BasisFunction Base;
        int ml;
    };
    
    

}} //namespace

#endif //_Atom_ml_Gaussian_BF_H_