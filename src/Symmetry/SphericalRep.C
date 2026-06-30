// File: Symmetry/SphericalRep.C  Representation of an orthogonal operation on a real-spherical shell.
//
// The spherical counterpart of CartesianShellRep.  A real solid harmonic chi_{l,m} is a fixed linear
// combination of the same-degree Cartesian monomials, so the (2l+1)x(2l+1) operation rep on the harmonic
// shell is obtained from the Cartesian rep by projecting through that c2s map.  The harmonic subspace is
// rotation-invariant, so the projection is EXACT (not a fit):
//     D_sph = (C C^T)^{-1} C  D_cart  C^T          ( C = the nSph x nCart coefficient matrix )
// with the same transform convention as CartesianShellRep: chi_m(R^{-1} r) = sum_n D(n,m) chi_n(r), and
// R |-> D_sph(R) a faithful representation.  Pure math: the c2s coefficients are supplied by the caller
// (the molecular basis side passes its SolidHarmonics), so qcSymmetry stays self-contained + LAPACK-free.
module;
#include <vector>
#include <array>
#include <utility>
export module qchem.Symmetry.SphericalRep;
export import qchem.Symmetry.CartesianRep;   // IVec3, CartesianShellRep (reused), rmat_t / Matrix3D

export namespace qchem::Symmetry
{

//! One harmonic's Cartesian expansion: a list of (monomial exponents, coefficient) terms.  c2s[m] is the
//! m-th real solid harmonic of the shell, in the BASIS'S OWN m-ordering (the convention must match the
//! basis these reps will adapt -- see the per-basis extractors).  Overall per-harmonic scale is irrelevant
//! (it cancels in the projection), so any convenient normalisation of the coefficients works.
using HarmonicC2S = std::vector<std::vector<std::pair<IVec3, double>>>;

//! Operation rep D_sph(R) on a real-spherical shell whose harmonics have the Cartesian expansion \a c2s.
rmat_t SphericalShellRep(const Matrix3D<double>& R, const HarmonicC2S& c2s);

} //namespace
