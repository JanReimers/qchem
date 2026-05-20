// File: SCFIterator/Types.C  Define types used throughout the SCFIterator module.
module;

export module qchem.SCFIterator.Types;
export import qchem.BasisSet;

export namespace qchem::SCFIterator
{
    using bs_t=BasisSet::BasisSet<double>;
}

