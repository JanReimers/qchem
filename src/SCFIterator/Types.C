// File: SCFIterator/Types.C  Define types used throughout the SCFIterator module.
module;

export module qchem.SCFIterator.Types;
export import qchem.BasisSet;

export namespace qchem::SCFIterator
{
    template <class T> using tbs_t=BasisSet::tBasisSet<T>;
    using rbs_t=tbs_t<double>;
}

