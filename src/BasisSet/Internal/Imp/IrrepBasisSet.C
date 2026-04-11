// File: TBasisSetImplementation.C
module;
// #include <tuple>
// #include <iostream>
// #include <vector>
// #include <memory>
#include <cassert>


module qchem.BasisSet.Internal.IrrepBasisSet;


template <class T> const Symmetry& IrrepBasisSet_Common<T>::GetSymmetry() const
{
    assert(itsSymmetry);
    return *itsSymmetry;
}

template <class T> Irrep_QNs IrrepBasisSet_Common<T>::GetIrrep(const Spin& s) const
{
    assert(itsSymmetry);
    return Irrep_QNs(s,itsSymmetry);
}


//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> Orbital_IBS_Common<T>::Orbital_IBS_Common()
{
};




template class IrrepBasisSet_Common<double>;
template class Orbital_IBS_Common<double>;

