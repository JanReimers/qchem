// File: Symmetry/CartesianRep.C  Representation of an orthogonal operation on a Cartesian shell.
//
// Stage 2 of the molecular-symmetry plan: how a Cartesian-Gaussian angular shell transforms
// under a point-group operation.  A Cartesian basis function has an angular monomial
// p(x,y,z) = x^nx y^ny z^nz of total degree L.  Under the operation with matrix R (about the
// centroid) a function transforms as phi(r) -> phi(R^{-1} r), so the monomial maps to
// p(R^{-1} u), a degree-L polynomial = a linear combination of the shell's monomials.  This
// builds that (#components x #components) matrix.  Pure math (no basis/Gaussian knowledge);
// the molecular basis rep-builder (in BasisSet1) combines it with the center permutation.
module;
#include <vector>
#include <array>
export module qchem.Symmetry.CartesianRep;
export import qchem.Types;       // rmat_t
export import qchem.Matrix3D;    // Matrix3D

export namespace Symmetry
{

using IVec3 = std::array<int,3>;   // Cartesian monomial exponents (nx,ny,nz)

// Representation matrix of the orthogonal operation R on a complete Cartesian shell whose
// components are the monomials `exps` (all of the same total degree L), in the given order.
// D(b,a) is the coefficient of monomial exps[b] in  p_a(R^{-1} u),  so the basis functions
// transform as  phi_a(R^{-1} r) = sum_b D(b,a) phi_b(r).  R |-> D(R) is a faithful
// representation: D(R1) D(R2) = D(R1 R2).  For L = 1 (a p-shell), D(R) = R.  Unnormalized:
// the molecular rep-builder applies the per-component Gaussian normalization separately.
rmat_t CartesianShellRep(const Matrix3D<double>& R, const std::vector<IVec3>& exps);

} //namespace
