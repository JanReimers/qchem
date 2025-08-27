// File: BSpline/IE.C Common IE code for BSpline basis sets.
module;
#include <bspline/Core.h>
#include <memory>
export module qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.Basisset.Atom.radial.BSpline.BFGrouper;
import qchem.Basisset.Atom.radial.BSpline.IEC;

import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.Cache4;
export import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.IEClient;

import qchem.BasisSet.Internal.ERI4;

export namespace BSpline
{

template <class T, size_t K> class IE_BS_2E 
    : public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper<K>
{
public:
    IE_BS_2E(BS_Evaluator* bse,AtomIE_BS_2E_Angular* a) : itsAngular(a), itsEvaluator(bse) {};
    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const ::IrrepIEClient*, IBS_Evaluator*);
private: 
    std::unique_ptr<AtomIE_BS_2E_Angular> itsAngular;
private: 
    BS_Evaluator* itsEvaluator;
};


} //namespace
