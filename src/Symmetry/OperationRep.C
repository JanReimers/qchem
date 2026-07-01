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
export module qchem.Symmetry.OperationRep;
export import qchem.Symmetry.CartesianRep;   // IVec3, CartesianShellRep, rmat_t / Matrix3D
export import qchem.Symmetry.SphericalRep;   // HarmonicC2S, SphericalShellRep

export namespace qchem::Symmetry
{

//---------------------------------------------------------------------------------------
// One angular shell as the representation builder sees it.  Shells related by a symmetry operation share
// the same shellType (same radial exponents + L), so the center permutation can match a shell to its
// image; offset is the shell's first index in the global AO ordering.  A shell is EITHER Cartesian (its
// components are the `monomials`) OR spherical (its components are the real solid harmonics whose
// Cartesian expansion is `c2s`, in the basis's own m-ordering) -- exactly one of the two is populated.
struct AoShell
{
    int                 shellType;
    rvec3_t             center;
    std::vector<IVec3>  monomials;   //!< Cartesian components, in AO order (empty for a spherical shell)
    std::vector<double> norm;        //!< per-component normalization N_a (size = #components)
    size_t              offset;
    HarmonicC2S         c2s = {};    //!< spherical harmonics' Cartesian expansion (empty for a Cartesian shell)

    bool   IsSpherical()  const { return !c2s.empty(); }
    size_t nComponents()  const { return IsSpherical() ? c2s.size() : monomials.size(); }
};

// Representation matrix of the operation R (acting about `origin`) on the whole AO basis described by
// `shells`: a permutation of shells combined with each shell's angular rep (Cartesian or spherical) and
// per-component normalization.  Entry M(I,J) is the coefficient of normalized AO I in the transform of
// normalized AO J, so phi_J(R^{-1} r) = sum_I M(I,J) phi_I(r).  R |-> M(R) is a representation:
// M(R1) M(R2) = M(R1 R2).
rmat_t BuildOperationRep(const std::vector<AoShell>& shells, const Matrix3D<double>& R,
                         const rvec3_t& origin, double tol);

} //namespace
