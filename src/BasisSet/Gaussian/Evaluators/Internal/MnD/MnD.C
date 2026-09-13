// File: BasisSet/Gaussian/Evaluators/Internal/MnD/MnD.C  Umbrella for the generic McMurchie-Davidson core.
//
// Convenience aggregator: a single `import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD` pulls in the
// whole generic MnD surface (Index3, Triangle3D, RNLM, and the Hermite recursion as it lands here), so
// client code -- PG_Cart_MnD today, PG_Spherical_MnD tomorrow -- needs only one import.
//
// This umbrella is itself `.Internal.`, so it is never re-exported across a library boundary (CLAUDE.md):
// the generic MnD core stays an internal implementation detail of the molecular basis-set library.
export module qchem.BasisSet.Gaussian.Evaluators.Internal.MnD;
export import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Index3;
export import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Triangle3D;
export import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.RNLM;

namespace qchem {

} // namespace qchem