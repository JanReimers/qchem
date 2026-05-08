// File: SCFIterator/Types.C  Define types used throughout the SCFIterator module.
module;

export module qchem.SCFIterator.Types;
#ifdef LegacyBasisSet

export import qchem.BasisSet;

export namespace qchem::SCFIterator
{
    using bs_t=BasisSet;
}
#else
export import qchem.BasisSet1;

export namespace qchem::SCFIterator
{
    using bs_t=BasisSet1::BasisSet<double>;
}

#endif