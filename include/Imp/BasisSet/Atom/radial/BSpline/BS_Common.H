// File: Atom/radial/BSpline/BS_Common.H  l/ml/kappa/mj independent part of BasisSet for atom BSpline Basis Sets.
#ifndef _Atom_radial_BSpline_BS_Common_H_
#define _Atom_radial_BSpline_BS_Common_H_

#include "../../../BS_Common.H"
#include "IE_HF.H"

namespace BSpline
{

// Common base handles all the radial aspects.
template <size_t K> class BS_Common
: public ::BS_Common
, public IE_BS_2E<double,K>
{
protected:
    virtual void Insert(bs_t* bs);
private:
    using BSpline::BFGrouper<K>::unique_spv;
    using BSpline::BFGrouper<K>::LMax;
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const;
};

}

#endif //_Atom_radial_BSpline_BS_Common_H_