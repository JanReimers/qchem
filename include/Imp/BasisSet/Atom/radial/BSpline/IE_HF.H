// File: BSpline/IE_HF.H  Integrals for BSpline HF basis sets.
#ifndef _BSpline_IE_HF_H_
#define _BSpline_IE_HF_H_

#include "IE.H"
#include "IEC.H"
#include "BFGrouper.H"
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet;
import qchem.BasisSet.Imp.HeapDB;
import qchem.BasisSet.ERI4;
import qchem.HF_IBS;

namespace BSpline
{

template <class T, size_t K> class IE_BS_2E 
    : public virtual ::AtomIE_BS_2E_Angular
    , public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper<K>
{
    typedef typename ::AtomIE_BS_2E_Angular::RVec RVec;
public:
    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const ::IrrepIEClient*);
};

} //namespace

#endif //_BSpline_IE_HF_H_