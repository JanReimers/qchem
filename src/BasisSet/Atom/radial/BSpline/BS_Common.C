// File: Atom/radial/BSpline/BS_Common.C  l/ml/kappa/mj independent part of BasisSet for atom BSpline Basis Sets.
module;
#include <vector>
#include <iostream>
export module qchem.Basisset.Atom.radial.BSpline.BS_Common; 
import qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.Basisset.Atom.radial.BSpline.BFGrouper;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet.Internal.IBS_Common;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.Cache4;
import oml.Vector;
export import qchem.Types;

export namespace BSpline
{
template <size_t K> class IrrepBasisSet
    : public virtual ::IrrepBasisSet
    , public         TIBS_Common1<double>
    , public         IrrepIEClient<K>
{
public:
    IrrepBasisSet(size_t Ngrid,double rmin, double rmax, Symmetry*,size_t L);
    IrrepBasisSet(size_t Ngrid,double rmin, double rmax, Symmetry*,size_t L, const std::vector<int>& ml);
    virtual std::ostream&  Write(std::ostream&    ) const;

};
template <size_t K> class BS_Common
: public ::BS_Common
, public IE_BS_2E<double,K>
{
protected:
    BS_Common() : itsRkCache(0) {};
    ~BS_Common();
    virtual void Insert(bs_t* bs);
    void BuildCache(size_t lmax);
private:
    using BSpline::BFGrouper<K>::unique_spv;
    using BSpline::BFGrouper<K>::LMax;
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const;

    RkCache<K>* itsRkCache;
};

}

