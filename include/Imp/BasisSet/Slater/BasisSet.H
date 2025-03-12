// File Slater/BasisSet.H
#ifndef _Slater_BasisSet_H_
#define _Slater_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"
#include "Imp/BasisSet/AtomIE.H"

namespace Slater
{

class BasisSet 
    : public ::BasisSetImp
    , public ::AtomIE_BS_2E<double>
{
protected:
    typedef ::AtomIE_BS_2E<double>::RVec RVec;
    typedef ::BasisSet Base;
    typedef Base::bs_t bs_t;
public:
    BasisSet() {};
    BasisSet(const LAParams&,size_t N, double minexp, double maxexp, size_t Lmax);
    virtual void Insert(bs_t* bs);

    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const;
};

} //namespace Slater

#endif // _Slater_BasisSet_H_
