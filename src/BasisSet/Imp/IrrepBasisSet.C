// File: BasisSet/Imp/IrrepBasisSet.C
module qchem.IrrepBasisSet;

template <class T> size_t IrrepBasisSet<T>::GetVectorSize() const 
{
    return GetNumFunctions();
}

// g++ 15.2 BUG very hard to auto instance this intermediary class.
template class IrrepBasisSet<double>;

