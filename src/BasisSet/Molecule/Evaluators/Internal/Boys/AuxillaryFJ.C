// File: BasisSet/Molecule/Evaluators/Internal/Boys/AuxillaryFJ.C  The Boys auxiliary function F_j(T).
//
// Not MnD-specific: every analytic Gaussian integral scheme (MnD, PRISM, OS, HGP, ...) needs the Boys
// function, so it lives in its own folder beside MnD.  Multiple implementations may land here later.
export module qchem.BasisSet.Molecule.Evaluators.Internal.Boys.AuxillaryFJ;
import qchem.Types;

export namespace qchem::BasisSet::Molecule::Evaluators::Internal::Boys
{
//
//  This class encapsulates the calculation of:
//
//          /1
//  Fj(T) = |  u**(2*j) * exp(-T*u**2)du
//          /0
//
class AuxillaryFJ
{
public:
    void GetFjAt(double T, rvec_t& Fj) const;
    static const int thejMax;

private:
    static double    theLookUp[121][23];
};

} //namespace

