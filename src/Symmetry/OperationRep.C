// File: Symmetry/OperationRep.C  Representation of an orthogonal operation on a whole AO basis.
//
// The angular-agnostic "whole-basis" layer above the per-shell reps: it combines the center
// permutation (a shell's center C -> R C, matched by shellType) with each shell's own angular rep --
// CartesianShellRep for a Cartesian-monomial shell, SphericalShellRep for a real-solid-harmonic shell --
// plus the per-component normalization.  Living above BOTH per-shell reps is what lets one AoShell list
// mix conventions and resolves the CartesianRep <-> SphericalRep dependency (SphericalRep reuses
// CartesianShellRep, so the dispatcher cannot live in either).
module;
#include <vector>
#include <memory>
export module qchem.Symmetry.OperationRep;
export import qchem.Symmetry.ShellRep;   // ShellRep (the angular rep abstraction), rmat_t / Matrix3D

export namespace qchem::Symmetry
{

//---------------------------------------------------------------------------------------
// One angular shell as the representation builder sees it: its place in the permutation structure (center,
// shellType so symmetry-equivalent shells match, offset in the global AO ordering) + per-component
// normalization + its angular rep.  The angular rep is a ShellRep (Cartesian or spherical -- the builder
// neither knows nor cares): AoShell owns no monomials/c2s and no angular-kind flag, so there is nothing to
// discriminate.  The per-basis extractor constructs the concrete ShellRep.
struct AoShell
{
    int                             shellType;
    rvec3_t                         center;
    std::vector<double>             norm;    //!< per-component normalization N_a (size = #components)
    size_t                          offset;
    std::shared_ptr<const ShellRep> rep;     //!< the shell's angular operation rep (Cartesian or spherical)

    size_t nComponents() const { return rep->nComponents(); }
};

// Representation matrix of the operation R (acting about `origin`) on the whole AO basis described by
// `shells`: a permutation of shells combined with each shell's angular rep (Cartesian or spherical) and
// per-component normalization.  Entry M(I,J) is the coefficient of normalized AO I in the transform of
// normalized AO J, so phi_J(R^{-1} r) = sum_I M(I,J) phi_I(r).  R |-> M(R) is a representation:
// M(R1) M(R2) = M(R1 R2).
rmat_t BuildOperationRep(const std::vector<AoShell>& shells, const Matrix3D<double>& R,
                         const rvec3_t& origin, double tol);

} //namespace
