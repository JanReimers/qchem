//! \file
//! \brief Abstract per-shell operation representation -- the angular capability an \c AoShell needs, as an
//! abstraction (Dependency Inversion): "given an orthogonal operation \f$R\f$, hand me the
//! \f$n_c\times n_c\f$ rep matrix on my components."  Concrete implementations know whether the components
//! are Cartesian monomials (\c CartesianShellRep) or real solid harmonics (\c SphericalShellRep); the
//! whole-basis rep-builder (\c BuildOperationRep) depends only on THIS interface, so it neither branches on
//! angular kind nor imports the concretes.
module;
export module qchem.Symmetry.Molecule.ShellRep;
export import qchem.Types;       // rmat_t, rmat3d_t (the fixed 3x3 operation matrix)

export namespace qchem::Symmetry::Molecule
{

//! \brief One shell's operation rep.  \c Rep(R) is the \f$n_c\times n_c\f$ matrix \f$D(b,a)\f$ with
//! \f$\phi_a(R^{-1}r)=\sum_b D(b,a)\,\phi_b(r)\f$ on the shell's components (in the shell's own order).
//! \f$R\mapsto D(R)\f$ is a faithful representation: \f$D(R_1)D(R_2)=D(R_1R_2)\f$.  Unnormalized -- the
//! whole-basis builder applies the per-component normalization separately.
class ShellRep
{
public:
    virtual ~ShellRep() {}
    virtual size_t nComponents() const = 0;
    virtual rmat_t Rep(const rmat3d_t& R) const = 0;
};

} //namespace
