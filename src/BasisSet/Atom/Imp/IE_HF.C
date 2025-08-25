// File: AtomIE.C Common HF IE code for all atom basis sets.
module;
#include <iostream>
#include <memory>
#include <cassert>
module qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IEClient;

template <class T> void AtomIE_BS_2E<T>::Append(const IrrepIEClient* ciec, IBS_Evaluator* eval)
{
    assert(ciec);
    DB_BS_2E<T>::Append(ciec);
    assert(eval);
    itsEvaluator->Register(eval);
}
template <class T> ERI4 AtomIE_BS_2E<T>::MakeDirect  (const IrrepIEClient* _a, const IrrepIEClient* _c) const
{
    const IBS_Evaluator* a=dynamic_cast<const IBS_Evaluator*>(_a);
    const IBS_Evaluator* c=dynamic_cast<const IBS_Evaluator*>(_c);
    assert(a);
    assert(c);
    return itsEvaluator->Direct(a,c);
}
template <class T> ERI4 AtomIE_BS_2E<T>::MakeExchange(const IrrepIEClient* _a, const IrrepIEClient* _c) const
{
    const IBS_Evaluator* a=dynamic_cast<const IBS_Evaluator*>(_a);
    const IBS_Evaluator* c=dynamic_cast<const IBS_Evaluator*>(_c);
    assert(a);
    assert(c);
    return itsEvaluator->Exchange(a,c);
}

template class AtomIE_BS_2E<double>;

 