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
template <size_t K> class BS_Common
: public ::BS_Common
, public IE_BS_2E<double,K>
{
protected:
    BS_Common(BS_Evaluator* eval, AtomIE_BS_2E_Angular* a) : IE_BS_2E<double,K>(eval,a), itsRkCache(0), itsEval(eval) {};
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
    BS_Evaluator* itsEval;
};

}

