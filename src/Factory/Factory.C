// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
export module qchem.Factory;

export import qchem.BasisSet.Radial.Factory;
export import qchem.BasisSet.Gaussian.Point.Factory;
export import qchem.BasisSet;

namespace qchem {
export using Real_BS=BasisSet::tBasisSet<double>;
export using Real_OIBS=BasisSet::Real_OIBS;
export using abs_t=BasisSet::Radial::Type;

} // namespace qchem