// File: Atom/radial/BSpline/BS_Common.C  l/ml/kappa/mj independent part of BasisSet for atom BSpline Basis Sets.
module;
#include <vector>
#include <iostream>
export module qchem.Basisset.Atom.radial.BSpline.BS_Common; 
import qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.BasisSet.Atom.Internal.radial.BSpline.Rk;
import qchem.Basisset.Atom.radial.BSpline.BFGrouper;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.IE;

import oml.Vector;
export import qchem.Types;

export namespace BSpline
{
template <size_t K> class IrrepBasisSet
    : public virtual Real_IBS
    , public         IrrepBasisSet_Common<double>
    , public         IrrepIEClient<K>
{
    typedef typename VectorFunction<double>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<double>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    typedef typename BSpline::IrrepIEClient<K> IEC;
    using IEC::splines;
    using IEC::rmin;
    using IEC::rmax;
    using AtomIrrepIEClient::ns;
public:
    IrrepBasisSet(size_t Ngrid,double rmin, double rmax, Symmetry*,size_t L);
    IrrepBasisSet(size_t Ngrid,double rmin, double rmax, Symmetry*,size_t L, const std::vector<int>& ml);
    virtual size_t  GetNumFunctions() const {return size();}
    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
    virtual std::ostream&  Write(std::ostream&    ) const;
};
template <size_t K> class BS_Common
: public ::BS_Common
, public IE_BS_2E<double,K>
{
protected:
    BS_Common(AtomIE_BS_2E_Angular* a) : IE_BS_2E<double,K>(a), itsRkCache(0) {};
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

