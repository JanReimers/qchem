// File: BasisSet/Molecule/Evaluators/PG_Spherical_MnD/SolidHarmonics.C
//
// The real-solid-harmonic -> Cartesian-monomial transform: the ONE thing a spherical-Gaussian basis
// needs beyond the Cartesian (PG_Cart_MnD) machinery.  A spherical shell function chi_{l,m} is a fixed
// linear combination of the same-L Cartesian monomials x^a y^b z^c (a+b+c = l), so every spherical
// integral is the corresponding Cartesian integral (GaussianRF kernels, already oracle-verified) summed
// with the coefficients returned here.
//
// SphericalShell(l) returns the 2l+1 functions ordered by m = -l..+l, each as its list of
// (Cartesian Polarization, coefficient) terms.  These are RAW monomial coefficients (degree-l harmonic
// polynomials); the evaluator folds in the per-component Cartesian normalisations and the overall
// spherical self-normalisation.  Within each m only the RELATIVE coefficients matter (the overall scale
// is fixed by normalising chi_{l,m} to unit self-overlap), so a convenient scale is used here.
module;
#include <vector>
export module qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD.SolidHarmonics;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;

export namespace BasisSet::Molecule::Evaluators::PG_Spherical_MnD
{
namespace Cart = ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;

// One term of a spherical function's Cartesian-monomial expansion.
struct CartTerm
{
    Cart::Polarization p;   // the Cartesian monomial x^p.n y^p.l z^p.m
    double             c;   // its coefficient
};

// The 2l+1 real solid harmonics for shell l, ordered m = -l..+l, each as its Cartesian expansion.
std::vector<std::vector<CartTerm>> SphericalShell(int l);

} //namespace
