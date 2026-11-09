// File: BasisSet/Gaussian/Evaluators/Internal/Boys/AuxillaryFJ.C  The Boys auxiliary function F_j(T).
//
// Not MnD-specific: every analytic Gaussian integral scheme (MnD, PRISM, OS, HGP, ...) needs the Boys
// function, so it lives in its own folder beside MnD.  Multiple implementations may land here later.
export module qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ;
import qchem.Types;

export namespace qchem::BasisSet::Gaussian::Evaluators::Internal::Boys
{
//
//  This class encapsulates the calculation of:
//
//          /1
//  Fj(T) = |  u**(2*j) * exp(-T*u**2)du
//          /0
//
//  The downward recursion starts from a 6-term Taylor expansion of the tabulated F_jmax, so the fixed table
//  (121 T points x 23 j columns) serves jmax <= thejMax = 16.  The bound is checked ONCE, at construction, and
//  an over-range jmax THROWS (CLAUDE.md "Where to check, and with what"): the old assert(Fj.size()<thejMax+6)
//  was both dead in Release and looser than the table (it let jmax reach 20, i.e. column 26 of 23).
//
class AuxillaryFJ
{
public:
    //! Throws std::out_of_range when jMax is outside [0, thejMax]: a caller with more angular momentum than the
    //! table holds (Ltotal = La+Lb+Lc+Ld > 16, e.g. four h shells) must be refused, not read past the table.
    explicit AuxillaryFJ(int jMax);
    //! Fills Fj[0..jMax]; Fj must have jMax+1 elements (the size the constructor checked).
    void GetFjAt(double T, rvec_t& Fj) const;
    static const int thejMax;

private:
    int              itsjMax;
    static double    theLookUp[121][23];
};

} //namespace

