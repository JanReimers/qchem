// File: Atom/ml/Slater_BF.H  r^l exp(-ar)*Y_lm type basis function 
#ifndef _Atom_ml_Slater_BF_H_
#define _Atom_ml_Slater_BF_H_

#include "../l/Slater_BF.H"

namespace Atom_ml
{
namespace Slater
{
    
class BasisFunction : public Atoml::Slater::BasisFunction
{
public:
    BasisFunction(double ex, int n, int l, int ml, double norm);
    
    virtual bool   operator==(const ::BasisFunction&) const;

    virtual BasisFunction* Clone(        ) const;
private:
    typedef Atoml::Slater::BasisFunction Base;
    int ml;
};


}} //namespace _Atom_ml_Slater_BF_H_

#endif //_SlaterBF_H_
