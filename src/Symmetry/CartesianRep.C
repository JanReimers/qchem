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

//---------------------------------------------------------------------------------------
// One Cartesian-Gaussian shell as the representation builder sees it.  Shells related by a
// symmetry operation share the same shellType (same radial exponents + L), so the center
// permutation can match a shell to its image; offset is the shell's first index in the
// global AO ordering.
struct AoShell
{
    int                 shellType;
    rvec3_t             center;
    std::vector<IVec3>  monomials;   // the shell's Cartesian components, in AO order
    std::vector<double> norm;        // per-component normalization N_a (same size as monomials)
    size_t              offset;
};

// Representation matrix of the operation R (acting about `origin`) on the whole AO basis
// described by `shells`: a permutation of shells (center C -> R C, matched by shellType)
// combined with the per-shell CartesianShellRep and per-component normalization.  Entry
// M(I,J) is the coefficient of normalized AO I in the transform of normalized AO J, so a
// normalized basis function transforms as phi_J(R^{-1} r) = sum_I M(I,J) phi_I(r).  As with
// the shell rep, R |-> M(R) is a representation: M(R1) M(R2) = M(R1 R2).
rmat_t BuildOperationRep(const std::vector<AoShell>& shells, const Matrix3D<double>& R,
                         const rvec3_t& origin, double tol);

} //namespace
