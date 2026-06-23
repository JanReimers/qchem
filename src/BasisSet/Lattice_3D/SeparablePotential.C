// File: BasisSet/Lattice_3D/SeparablePotential.C  Separable (Kleinman-Bylander) NONLOCAL potential.
//
// Rung 2 of "lineage A" (see doc/PlaneWavePlan.md): a norm-conserving pseudopotential's nonlocal part
// in Kleinman-Bylander separable form,
//     V_NL = Sum_{atom a} Sum_{projector p} |beta^a_p> D_p <beta^a_p| .
// In a plane-wave basis each projector contributes a reciprocal-space radial form factor
// beta-tilde_p(|q|) (q = k+G) and a KB coefficient D_p; the basis set (PlaneWave_IBS::
// MakeSeparablePotential) folds in 1/Omega and the structure-factor phase e^{-i(G-G').tau_a}.
//
// This is the nonlocal sibling of LocalPotential: the external one-body potential is V = V_loc +
// V_nonlocal, both model-parameterized contributions to the SAME external block.
//
// Angular channels: each projector carries an angular momentum l (AngularMomentum); MakeSeparablePotential
// weights its rank-1 radial product by the addition-theorem factor (2l+1) P_l(cos gamma), gamma the angle
// between k+G and k+G' -- the SAME angular structure the APW/LAPW sphere terms use.  l=0 (P_0=1) is the
// spherically-symmetric s-channel.  The radial form factor beta-tilde_l(|q|) is a model input here (the
// production norm-conserving / PAW route obtains it from a spherical-Bessel transform j_l of beta_l(r)).
module;
#include <cmath>

export module qchem.BasisSet.Lattice_3D.SeparablePotential;
import qchem.Math;

export namespace BasisSet::Lattice_3D
{

//! \brief A separable Kleinman-Bylander nonlocal potential, supplying per-species projector form
//! factors and KB coefficients.  Open/closed extension point for the nonlocal part of a pseudopotential.
class SeparablePotential
{
public:
    virtual ~SeparablePotential() {}
    //! Number of (radial) projectors for nuclear species \a Z.
    virtual size_t NumProjectors(int Z) const=0;
    //! Kleinman-Bylander coefficient \f$D_p\f$ for projector \a p of species \a Z.  [energy]
    virtual double Coefficient  (int Z, size_t p) const=0;
    //! Reciprocal-space radial projector form factor \f$\tilde\beta_p(|q|)\f$, \f$q=|k+G|\f$.
    virtual double Projector    (int Z, size_t p, double q) const=0;
    //! Angular momentum \a l of projector \a p (default 0 = s-channel).  Sets the \f$(2l+1)P_l(\cos\gamma)\f$
    //! angular weight of this projector's contribution.
    virtual int    AngularMomentum(int /*Z*/, size_t /*p*/) const {return 0;}
};

//! \brief A single Gaussian Kleinman-Bylander projector in channel \a l,
//! \f$\tilde\beta(q)=e^{-\sigma^2 q^2/2}\f$ with coefficient \f$D\f$ -- a minimal analytic demonstrator
//! of the separable nonlocal structure (each atom contributes a rank-1 \f$V_{NL}\f$ per channel).
//! Defaults to l=0 (the s-channel), so existing callers are unaffected.
class GaussianProjector : public SeparablePotential
{
public:
    GaussianProjector(double sigma, double D, int l=0) : itsSigma(sigma), itsD(D), itsL(l) {}
    virtual size_t NumProjectors(int) const {return 1;}
    virtual double Coefficient  (int, size_t) const {return itsD;}
    virtual double Projector    (int, size_t, double q) const {return std::exp(-0.5*itsSigma*itsSigma*q*q);}
    virtual int    AngularMomentum(int, size_t) const {return itsL;}
private:
    double itsSigma; //!< Projector width (Bohr).
    double itsD;     //!< KB coefficient (energy).
    int    itsL;     //!< Angular-momentum channel.
};

} //namespace
