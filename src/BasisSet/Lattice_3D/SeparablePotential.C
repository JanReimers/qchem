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
// Demonstrator scope: s-channel (l=0) projectors -- the form factor is the spherically-symmetric radial
// transform, with no angular dependence.  l>0 additionally needs Y_lm(q-hat) and a spherical-Bessel
// transform j_l (the production norm-conserving / PAW extension); the interface leaves room for it.
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
    //! Reciprocal-space radial projector form factor \f$\tilde\beta_p(|q|)\f$, \f$q=k+G\f$ (l=0).
    virtual double Projector    (int Z, size_t p, double q) const=0;
};

//! \brief A single s-channel (l=0) Gaussian Kleinman-Bylander projector,
//! \f$\tilde\beta(q)=e^{-\sigma^2 q^2/2}\f$ with coefficient \f$D\f$ -- a minimal analytic demonstrator
//! of the separable nonlocal structure (each atom contributes a rank-1 \f$V_{NL}\f$).
class GaussianProjector : public SeparablePotential
{
public:
    GaussianProjector(double sigma, double D) : itsSigma(sigma), itsD(D) {}
    virtual size_t NumProjectors(int) const {return 1;}
    virtual double Coefficient  (int, size_t) const {return itsD;}
    virtual double Projector    (int, size_t, double q) const {return std::exp(-0.5*itsSigma*itsSigma*q*q);}
private:
    double itsSigma; //!< Projector width (Bohr).
    double itsD;     //!< KB coefficient (energy).
};

} //namespace
