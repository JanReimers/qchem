// File: BasisSet/Molecule/Evaluators/Internal/MnD/Index3.C
//
// Generic 3-index (N,L,M) for the McMurchie-Davidson Hermite machinery.  This is the seam that keeps the
// generic MnD core (RNLM, the Hermite recursion, the Boys function) independent of any *angular*
// representation: RNLM and friends are addressed by a plain (N,L,M) Hermite index, NOT by the Cartesian
// PolarizedGaussian::Polarization.  The Cartesian layer converts: Polarization IS implicitly an Index3
// (it adds the monomial x^n y^l z^m on top), and a future Spherical layer would supply its own indices.
export module qchem.BasisSet.Molecule.Evaluators.Internal.MnD.Index3;

export namespace BasisSet::Molecule::Evaluators::Internal::MnD
{

struct Index3
{
    int n=0, l=0, m=0;
};

} //namespace
