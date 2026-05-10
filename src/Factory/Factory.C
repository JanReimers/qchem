// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
export module qchem.Factory;

#ifdef LegacyBasisSet
export import qchem.BasisSet;
export import qchem.BasisSet.Atom.Factory;
export import qchem.BasisSet.Molecule.Factory;
export namespace BasisSetAtomFactory=BasisSetAtom;
export using bs_t=BasisSet;
#else
export import qchem.BasisSet1.Atom.Factory;
export import qchem.BasisSet1.Molecule.Factory;
export import qchem.BasisSet1;
export using bs_t=BasisSet1::BasisSet<double>;
export using BasisSet=BasisSet1::BasisSet<double>;
export using Real_OIBS=BasisSet1::Real_OIBS;
export namespace BasisSetAtomFactory=BasisSet1::Atom;
#endif
