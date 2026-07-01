// File: Symmetry/ShellRep.C  Abstract per-shell operation representation.
//
// The angular capability an AoShell needs, as an abstraction (Dependency Inversion): "given an orthogonal
// operation R, hand me the (nc x nc) rep matrix on my components."  Concrete implementations know whether
// the components are Cartesian monomials (CartesianShellRep) or real solid harmonics (SphericalShellRep);
// the whole-basis rep-builder (BuildOperationRep) depends only on THIS interface, so it neither branches on
// angular kind nor imports the concretes.
module;
export module qchem.Symmetry.ShellRep;
export import qchem.Types;       // rmat_t
export import qchem.Matrix3D;    // Matrix3D

export namespace qchem::Symmetry
{

//! One shell's operation rep.  Rep(R) is the (nComponents x nComponents) matrix D(b,a) with
//! phi_a(R^{-1} r) = sum_b D(b,a) phi_b(r) on the shell's components (in the shell's own order).
//! R |-> D(R) is a faithful representation: D(R1) D(R2) = D(R1 R2).  Unnormalized -- the whole-basis
//! builder applies the per-component normalization separately.
class ShellRep
{
public:
    virtual ~ShellRep() {}
    virtual size_t nComponents() const = 0;
    virtual rmat_t Rep(const Matrix3D<double>& R) const = 0;
};

} //namespace
