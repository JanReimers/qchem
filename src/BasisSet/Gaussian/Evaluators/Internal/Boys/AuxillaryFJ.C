// File: BasisSet/Gaussian/Evaluators/Internal/Boys/AuxillaryFJ.C  The Boys auxiliary function F_j(T).
//
// Not MnD-specific: every analytic Gaussian integral scheme (MnD, PRISM, OS, HGP, ...) needs the Boys
// function, so it lives in its own folder beside MnD.  Multiple implementations may land here later.
module;
#include <cassert>
#include <cmath>
export module qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ;
import qchem.Types;
import qchem.CMath;    // Pi
import qchem.IntPow;   // uintpow

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
    void GetFjAt(double T, rvec_t& Fj) const {GetFjAt<rvec_t>(T,Fj);}
    //! The same into any indexable \a V.  INLINE (a module's out-of-line body is not importable): the cache-free nuclear kernel Enzyme
    //! differentiates (GaussianRF::NuclearGradient) evaluates F_j through it.  ⚠ uintpow's x == 0 guard makes the Taylor step's dt^1 term
    //! read a zero derivative exactly ON a table node; the only node a geometry reaches exactly is T = 0, where dT/dR = 0 too.
    template <class V> inline void GetFjAt(double T, V& Fj) const;
    static const int thejMax;
    static const double theLookUp[121][23];   //!< const: an Enzyme-differentiated read sees an inactive global

private:
    int              itsjMax;
};

template <class V> inline void AuxillaryFJ::GetFjAt(double T, V& Fj) const
{
    constexpr double Factorial[]={1,1,2,6,24,120,720};
    assert(size_t(Fj.size())==size_t(itsjMax+1)); // the size the constructor checked
    assert(T>=0                          );
    int jmax=itsjMax;
    double et=std::exp(-T);
    if(T<12)
    {
        int it=(int)(T*10.0+0.5);
        double Tstar=(double)it*0.1;
        double dt=Tstar-T;
        const double* FT=&theLookUp[it][0];
        double fjmax=FT[jmax];
        for (int k=1; k<=6; k++) fjmax+=FT[k+jmax]*uintpow(dt,k)/Factorial[k];
        Fj[jmax]=fjmax;
        for (int j=jmax-1; j>=0; j--) Fj[j]=(2*T*Fj[j+1]+et)/(2*j+1);
    }
    if(T>=12 && T < 30.0)
    {
        double g=0;
        if(         T<15) g = 0.4999489092 - 0.2473631686/T + 0.321180909/(T*T) - 0.3811559346/(T*T*T);
        if(T>=15 && T<18) g = 0.4998436875 - 0.24249438  /T + 0.24642845 /(T*T);
        if(T>=18 && T<24) g = 0.499093162  - 0.2152832   /T;
        if(T>=24        ) g = 0.49;
        assert(g!=0.0); // g never got assigned.
        Fj[0] = 0.5*std::sqrt(Pi/T) - et*g/T;
        for (int j=0; j<jmax; j++) Fj[j+1]=((2*j+1)*Fj[j]-et)/(2*T);
    }
    if (T>=30.0)
    {
        Fj[0]=0.5*std::sqrt(Pi/T);
        for (int j=0; j<jmax; j++) Fj[j+1]=((2*j+1)*Fj[j]-et)/(2*T);
    }
}

} //namespace

