// File: Symmetry/CartesianRep.C  Representation of an orthogonal operation on a Cartesian shell.
//
// Stage 2 of the molecular-symmetry plan: how a Cartesian-Gaussian angular shell transforms
// under a point-group operation.  A Cartesian basis function has an angular monomial
// p(x,y,z) = x^nx y^ny z^nz of total degree L.  Under the operation with matrix R (about the
// centroid) a function transforms as phi(r) -> phi(R^{-1} r), so the monomial maps to
// p(R^{-1} u), a degree-L polynomial = a linear combination of the shell's monomials.  This
// builds that (#components x #components) matrix.  Pure math (no basis/Gaussian knowledge);
// the molecular basis rep-builder (in BasisSet) combines it with the center permutation.
module;
#include <vector>
#include <array>
#include <utility>
export module qchem.Symmetry.CartesianRep;
export import qchem.Symmetry.ShellRep;   // ShellRep (the abstraction this implements), rmat_t / Matrix3D

export namespace qchem::Symmetry
{

using IVec3 = std::array<int,3>;   // Cartesian monomial exponents (nx,ny,nz)

//! The operation rep of a complete Cartesian shell whose components are the monomials `exps` (all of the
//! same total degree L), in the given order.  Rep(R)(b,a) is the coefficient of monomial exps[b] in
//! p_a(R^{-1} u), so phi_a(R^{-1} r) = sum_b Rep(b,a) phi_b(r).  R |-> Rep(R) is faithful: for L = 1 (a
//! p-shell) Rep(R) = R.  Unnormalized -- the whole-basis builder applies per-component normalization.
class CartesianShellRep : public ShellRep
{
public:
    explicit CartesianShellRep(std::vector<IVec3> exps) : itsExps(std::move(exps)) {}
    virtual size_t nComponents() const {return itsExps.size();}
    virtual rmat_t Rep(const Matrix3D<double>& R) const;
private:
    std::vector<IVec3> itsExps;
};

} //namespace
