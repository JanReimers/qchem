// File SphericalGaussian/BasisSet.H
#ifndef _SphericalGaussian_BasisSet_H_
#define _SphericalGaussian_BasisSet_H_

#include "../TCommon.H"
#include "../Atom_IE_HF.H"

namespace SphericalGaussian
{

// Common base handles all the radial aspects.
class BS_Common
: public ::BasisSetImp
, public ::AtomIE_BS_2E<double>
{
protected:
    virtual void Insert(bs_t* bs);
private:
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const;
};

class BasisSet 
: public BS_Common
, public ::AtomIE_BS_2E_l //Pick angular integrals.
{
protected:
    typedef ::AtomIE_BS_2E_Angular::RVec RVec;
    typedef ::BasisSet Base;
    typedef Base::bs_t bs_t;
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax);
};

} //namespace

#endif // _SphericalGaussian_BasisSet_H_
