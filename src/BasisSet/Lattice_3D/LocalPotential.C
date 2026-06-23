// File: BasisSet/Lattice_3D/LocalPotential.C  One-body LOCAL external potentials for plane-wave work.
//
// A local external potential, in a plane-wave basis, is fully specified by its reciprocal-space radial
// form factor v(|G|^2) per nuclear species.  The plane-wave basis set folds in the cell volume 1/Omega,
// the structure factor Sum_a e^{-iG.tau_a}, and the G=0 handling (see PlaneWave_IBS::MakeLocalPotential);
// a LocalPotential supplies ONLY the species form factor.
//
// This is the open/closed extension point for "lineage A" (plane waves + pseudopotentials, see
// doc/PlaneWavePlan.md): the bare Coulomb nucleus, the Gaussian-smeared nucleus (rung-1 local
// pseudopotential), and -- later -- tabulated norm-conserving pseudopotential form factors are all
// just implementations of FormFactor(); the assembler never changes.
module;
#include <cmath>

export module qchem.BasisSet.Lattice_3D.LocalPotential;
import qchem.Math; // FourPi, Pi

export namespace BasisSet::Lattice_3D
{

//! \brief A one-body local external potential, defined by its reciprocal-space form factor.
class LocalPotential
{
public:
    virtual ~LocalPotential() {}
    //! \brief \f$v(Z,|G|^2)\f$: the Fourier transform of a single atom's potential (\f$|G|>0\f$),
    //! excluding the \f$1/\Omega\f$ and the structure-factor phase (those are geometry, applied by
    //! the basis set).  [energy x volume]
    virtual double FormFactor(int Z, double G2) const=0;
};

//! \brief Bare nuclear Coulomb \f$v(G) = -4\pi Z/G^2\f$.  Physically exact but the 1s cusp makes the
//! plane-wave energy converge very slowly in \f$E_{cut}\f$.
class BareCoulomb : public LocalPotential
{
public:
    virtual double FormFactor(int Z, double G2) const {return -FourPi*Z/G2;}
};

//! \brief Gaussian-smeared nucleus (rung-1 "local pseudopotential"): the point charge is spread into a
//! Gaussian of width \f$\sigma\f$, giving \f$v(G) = -4\pi Z\,e^{-\sigma^2 G^2/2}/G^2\f$.  The smooth
//! charge removes the cusp, so the plane-wave energy converges rapidly with \f$E_{cut}\f$;
//! \f$\sigma\to 0\f$ recovers BareCoulomb.
class GaussianSmearedNucleus : public LocalPotential
{
public:
    explicit GaussianSmearedNucleus(double sigma) : itsSigma(sigma) {}
    virtual double FormFactor(int Z, double G2) const
    {
        return -FourPi*Z*std::exp(-0.5*itsSigma*itsSigma*G2)/G2;
    }
private:
    double itsSigma; //!< Smearing width (Bohr).
};

//! \brief Local part of a real norm-conserving pseudopotential in the analytic Goedecker / Hartwigsen-
//! Goedecker-Hutter (HGH) form [Hartwigsen, Goedecker, Hutter, PRB 58, 3641 (1998); Goedecker, Teter,
//! Hutter, PRB 54, 1703 (1996)].  In real space
//! \f$V_{loc}(r) = -\frac{Z_{ion}}{r}\,\mathrm{erf}\!\big(\tfrac{r}{\sqrt2\,r_{loc}}\big)
//!                + e^{-r^2/2r_{loc}^2}\sum_{i=1}^{4} C_i\,(r/r_{loc})^{2i-2}\f$:
//! the erf softens the \f$-Z_{ion}/r\f$ singularity (long-range Coulomb preserved, core pseudized),
//! and the Gaussian-polynomial fits the rest.  Both pieces Fourier-transform in closed form, so the
//! reciprocal form factor below is analytic -- no radial tables.  \f$Z_{ion}\f$ is the VALENCE charge.
//! Unlike the all-electron LAPW, absolute levels are shifted by the dropped \f$G=0\f$ term, but the
//! softness (fast \f$E_{cut}\f$ convergence) and band-energy DIFFERENCES are physical.
class HGH_LocalPotential : public LocalPotential
{
public:
    //! \a Zion = valence charge, \a rloc = local radius, \a c = {C1,C2,C3,C4} polynomial coefficients.
    HGH_LocalPotential(double Zion, double rloc, double c1, double c2, double c3=0.0, double c4=0.0)
        : itsZion(Zion), itsRloc(rloc), itsC1(c1), itsC2(c2), itsC3(c3), itsC4(c4) {}

    //! Hydrogen, HGH (LDA): no nonlocal projectors (H has no core), so this local part is the whole PP.
    static HGH_LocalPotential Hydrogen() {return HGH_LocalPotential(1.0, 0.2, -4.0663326, 0.6778322);}

    virtual double FormFactor(int /*Z*/, double G2) const     // Z ignored: itsZion is the species
    {
        double t=G2*itsRloc*itsRloc;                          // (G r_loc)^2
        double g=std::exp(-0.5*t);
        double coulomb=-FourPi*itsZion/G2 * g;                // softened -Z_ion/r tail
        double poly=itsC1 + itsC2*(3-t) + itsC3*(15-10*t+t*t) + itsC4*(105-105*t+21*t*t-t*t*t);
        double twopi32=std::pow(2*Pi, 1.5);
        return coulomb + twopi32*itsRloc*itsRloc*itsRloc * g * poly;
    }
private:
    double itsZion, itsRloc, itsC1, itsC2, itsC3, itsC4;
};

} //namespace
