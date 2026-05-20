// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
export module qchem.Factory;

export import qchem.BasisSet.Atom.Factory;
export import qchem.BasisSet.Molecule.Factory;
export import qchem.BasisSet;
export using bs_t=BasisSet::BasisSet<double>;
export using Real_BS=BasisSet::BasisSet<double>;
export using Real_OIBS=BasisSet::Real_OIBS;
export namespace BasisSetAtomFactory=BasisSet::Atom;
