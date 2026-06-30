// File: BasisSet/LocalPotential.C  One-body LOCAL external potentials (plane-wave / pseudopotential).
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
#include <functional>
#include <map>
#include <memory>
#include <cassert>

export module qchem.Pseudopotential.LocalPotential;
import qchem.Math; // FourPi, Pi

export namespace qchem::Pseudopotential
{

//! \brief The RECIPROCAL-space face of a local external potential: its per-species form factor.  A radial
//! function has two spectral views; this is the \f$\tilde v(G^2)\f$ one that a PLANE-WAVE basis consumes.
class LocalPotential_Q
{
public:
    virtual ~LocalPotential_Q() {}
    //! \brief \f$v(Z,|G|^2)\f$: the Fourier transform of a single atom's potential (\f$|G|>0\f$),
    //! excluding the \f$1/\Omega\f$ and the structure-factor phase (those are geometry, applied by
    //! the basis set).  [energy x volume]
    virtual double FormFactor(int Z, double G2) const=0;

    //! \brief The finite \f$G\to 0\f$ limit of \f$v(G)+4\pi Z/G^2\f$, i.e. \f$\alpha=\int[V_{loc}(r)+Z/r]\,d^3r\f$
    //! (the form factor with its divergent Coulomb tail removed).  This is the "\f$\alpha\f$" / G=0
    //! alignment constant a plane-wave total energy needs: the \f$G=0\f$ electron-ion energy is
    //! \f$(N_{el}/\Omega)\sum_a\alpha_a\f$, the uniform shift dropped along with the \f$G=0\f$ potential.
    //! Default 0 (a pure Coulomb tail has no finite remainder). [energy x volume]
    virtual double FormFactorG0(int Z) const {return 0.0;}
};

//! \brief The REAL-space face: the same radial function as \f$V_{loc}(Z,r)\f$, the view a MOLECULAR /
//! ATOMIC (Gaussian or radial) basis consumes by quadrature \f$\langle\chi_i|V_{loc}|\chi_j\rangle\f$.
//! Split from the reciprocal face (the AO/FT axis, one level down in the PP model -- see
//! doc/MolecularPseudopotentialPlan.md section 2): a reciprocal-only model is-a LocalPotential_Q and is
//! simply not passable where a real-space view is required -- the wall is a type, not an NA-assert.
class LocalPotential_R
{
public:
    virtual ~LocalPotential_R() {}
    //! \brief \f$V_{loc}(Z,r)\f$ in a.u. (the full one-atom local potential; finite at \f$r=0\f$ for a
    //! pseudopotential -- no nuclear cusp).
    virtual double Vloc(int Z, double r) const=0;
};

//! \brief A one-body local external potential with BOTH spectral views (the closed-form models have both)
//! plus its ion charge.  The diamond LocalPotential_Q + LocalPotential_R is harmless (the project's house
//! style); the basis selects the view matching its space.
class LocalPotential
    : public virtual LocalPotential_Q
    , public virtual LocalPotential_R
{
public:
    //! \brief The ION charge this local potential's \f$-Z_{ion}/r\f$ tail carries -- the charge the ion-ion
    //! (Ewald) sum needs.  Default = the (true) nuclear charge \a Z (all-electron: BareCoulomb, Gaussian-
    //! smeared).  A pseudopotential overrides with its VALENCE charge (HGH: \f$Z_{ion}\f$), Z-independent.
    //! This is why it is a callback, not a getter: the answer depends on the model, not just on Z.
    virtual double Zion(int Z) const {return double(Z);}

    // --- Adapt Zion to the plain CALLBACK the ion-ion (Ewald) term consumes.  The reciprocal-space
    //     assembly (Integrals_Pseudo) now takes the abstract model DIRECTLY, but Ewald (NuclearRepulsion) lives
    //     BELOW the pseudopotential layer in qcStructure, so it must take a neutral std::function. ---
    std::function<double(int)>        ZionFn()         const {return [this](int Z){return Zion(Z);};}
};

//! \brief Bare nuclear Coulomb \f$v(G) = -4\pi Z/G^2\f$.  Physically exact but the 1s cusp makes the
//! plane-wave energy converge very slowly in \f$E_{cut}\f$.
class BareCoulomb : public LocalPotential
{
public:
    virtual double FormFactor(int Z, double G2) const {return -FourPi*Z/G2;}
    virtual double Vloc      (int Z, double r)  const {return -double(Z)/r;}   // -Z/r (singular at r=0)
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
    //! \f$v(G)+4\pi Z/G^2 = 4\pi Z(1-e^{-\sigma^2G^2/2})/G^2 \to 2\pi Z\sigma^2\f$ as \f$G\to0\f$.
    virtual double FormFactorG0(int Z) const {return 2*Pi*Z*itsSigma*itsSigma;}
    //! Real space: the point charge smeared into a Gaussian -> \f$V(r)=-Z\,\mathrm{erf}(r/\sqrt2\sigma)/r\f$
    //! (finite \f$-Z\sqrt{2/\pi}/\sigma\f$ at \f$r=0\f$).
    virtual double Vloc(int Z, double r) const
    {
        double x=r/(std::sqrt(2.0)*itsSigma);
        return (r>1e-12) ? -double(Z)*std::erf(x)/r
                         : -double(Z)*std::sqrt(2.0/Pi)/itsSigma;   // erf(x)/x -> 2/sqrt(pi) as x->0
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
    //! Real per-element parameters come from the GTH database via GetGTH (GTH_Potentials.C), not
    //! hardcoded factories.
    HGH_LocalPotential(double Zion, double rloc, double c1, double c2, double c3=0.0, double c4=0.0)
        : itsZion(Zion), itsRloc(rloc), itsC1(c1), itsC2(c2), itsC3(c3), itsC4(c4) {}

    virtual double FormFactor(int /*Z*/, double G2) const     // Z ignored: itsZion is the species
    {
        double t=G2*itsRloc*itsRloc;                          // (G r_loc)^2
        double g=std::exp(-0.5*t);
        double coulomb=-FourPi*itsZion/G2 * g;                // softened -Z_ion/r tail
        double poly=itsC1 + itsC2*(3-t) + itsC3*(15-10*t+t*t) + itsC4*(105-105*t+21*t*t-t*t*t);
        double twopi32=std::pow(2*Pi, 1.5);
        return coulomb + twopi32*itsRloc*itsRloc*itsRloc * g * poly;
    }
    //! Real space (the closed-form inverse FT documented above):
    //! \f$V_{loc}(r)=-\tfrac{Z_{ion}}{r}\mathrm{erf}(\tfrac{r}{\sqrt2 r_{loc}})
    //! + e^{-r^2/2r_{loc}^2}\sum_i C_i (r/r_{loc})^{2i-2}\f$.  Finite at \f$r=0\f$ (no nuclear cusp).
    virtual double Vloc(int /*Z*/, double r) const
    {
        double x=r/itsRloc, x2=x*x;
        double gpoly=std::exp(-0.5*x2)*(itsC1 + itsC2*x2 + itsC3*x2*x2 + itsC4*x2*x2*x2);
        double a=std::sqrt(2.0)*itsRloc;                       // erf argument scale
        double erfterm = (r>1e-12) ? -itsZion*std::erf(r/a)/r
                                   : -itsZion*2.0/(std::sqrt(Pi)*a);   // erf(r/a)/r -> 2/(sqrt(pi) a)
        return erfterm + gpoly;
    }
    //! G=0 alignment \f$\alpha=\int[V_{loc}+Z_{ion}/r]\,d^3r = 2\pi Z_{ion}r_{loc}^2
    //! + (2\pi)^{3/2}r_{loc}^3(C_1+3C_2+15C_3+105C_4)\f$ (the \f$G\to0\f$ finite part: the softened
    //! Coulomb leaves \f$2\pi Z r_{loc}^2\f$, the Gaussian-polynomial leaves the moment sum).
    virtual double FormFactorG0(int /*Z*/) const
    {
        double twopi32=std::pow(2*Pi, 1.5);
        double moments=itsC1 + 3*itsC2 + 15*itsC3 + 105*itsC4;
        return 2*Pi*itsZion*itsRloc*itsRloc + twopi32*itsRloc*itsRloc*itsRloc*moments;
    }
    virtual double Zion(int /*Z*/) const override {return itsZion;}   // valence charge (Z-independent)
private:
    double itsZion, itsRloc, itsC1, itsC2, itsC3, itsC4;
};

//! \brief A multi-species local potential: a router keyed by atomic number \a Z that forwards each
//! query to the per-species sub-model.  This is ALL that multi-species (ionic) crystals need on the
//! local side -- every LocalPotential method already takes \a Z (single-species HGH ignores it; this
//! USES it to dispatch), so the basis assembly (which calls FormFactor(a->itsZ,...) per atom) is
//! unchanged.  Hand one of these to the external term and NaF / CsI just work.
class MultiSpecies_LocalPotential : public LocalPotential
{
public:
    //! Register species \a Z's local model (atomic number Z, e.g. 11 for Na -- the atoms' itsZ).
    void Add(int Z, std::shared_ptr<const LocalPotential> model) {itsByZ[Z]=std::move(model);}
    virtual double FormFactor  (int Z, double G2) const override {return Get(Z).FormFactor(Z,G2);}
    virtual double FormFactorG0(int Z)            const override {return Get(Z).FormFactorG0(Z);}
    virtual double Vloc        (int Z, double r)  const override {return Get(Z).Vloc(Z,r);}
    virtual double Zion        (int Z)            const override {return Get(Z).Zion(Z);}
private:
    const LocalPotential& Get(int Z) const
    {
        auto it=itsByZ.find(Z);
        assert(it!=itsByZ.end() && "MultiSpecies_LocalPotential: no model registered for this species Z");
        return *it->second;
    }
    std::map<int, std::shared_ptr<const LocalPotential>> itsByZ;   //!< atomic number -> that species' local model
};

} //namespace
