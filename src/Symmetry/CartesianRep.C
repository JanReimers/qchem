//! \file
//! \brief Representation of an orthogonal operation on a Cartesian shell -- how a Cartesian-Gaussian angular
//! shell transforms under a point-group operation.  A Cartesian basis function has an angular monomial
//! \f$p(x,y,z)=x^{n_x}y^{n_y}z^{n_z}\f$ of total degree \f$L\f$.  Under the operation with matrix \f$R\f$
//! (about the centroid) a function transforms as \f$\phi(r)\to\phi(R^{-1}r)\f$, so the monomial maps to
//! \f$p(R^{-1}u)\f$, a degree-\f$L\f$ polynomial = a linear combination of the shell's monomials.  This
//! builds that \f$n_c\times n_c\f$ matrix.  Pure math (no basis/Gaussian knowledge); the molecular basis
//! rep-builder (in BasisSet) combines it with the centre permutation.
module;
#include <vector>
#include <array>
#include <utility>
export module qchem.Symmetry.CartesianRep;
export import qchem.Symmetry.ShellRep;   // ShellRep (the abstraction this implements), rmat_t / Matrix3D

export namespace qchem::Symmetry
{

using IVec3 = std::array<int,3>;   //!< Cartesian monomial exponents \f$(n_x,n_y,n_z)\f$

//! \brief The operation rep of a complete Cartesian shell whose components are the monomials \a exps (all of
//! the same total degree \f$L\f$), in the given order.  \f$\text{Rep}(R)(b,a)\f$ is the coefficient of
//! monomial \c exps[b] in \f$p_a(R^{-1}u)\f$, so \f$\phi_a(R^{-1}r)=\sum_b \text{Rep}(b,a)\,\phi_b(r)\f$.
//! \f$R\mapsto\text{Rep}(R)\f$ is faithful: for \f$L=1\f$ (a p-shell) \f$\text{Rep}(R)=R\f$.  Unnormalized --
//! the whole-basis builder applies per-component normalization.
class CartesianShellRep : public ShellRep
{
public:
    explicit CartesianShellRep(std::vector<IVec3> exps) : itsExps(std::move(exps)) {}
    virtual size_t nComponents() const {return itsExps.size();}
    virtual rmat_t Rep(const rmat3d_t& R) const;
private:
    std::vector<IVec3> itsExps;
};

} //namespace
