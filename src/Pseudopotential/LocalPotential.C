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
#include <vector>
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
    // --- The local potential is the sum of a LONG-range (softened-Coulomb / Gaussian core-charge) part and a
    //     SHORT-range (poly x Gaussian) remainder -- the CP2K local-PP split (doc/GPWPlan.md 0e-PP).  This
    //     decomposition is PRIMARY: a model supplies its own \c FormFactorLong (+ \c FormFactorShort if it has
    //     a compact core), and the full \f$v(G)\f$ + G=0 alignment are their sums, provided here.  The split is
    //     what the periodic KS assembly consumes -- the LONG part folds into the G-space Poisson (the Hartree
    //     term, one electrostatics solve), the SHORT part is the external term's compact remainder. ---

    //! \brief LONG-range part of \f$v(Z,|G|^2)\f$ (\f$|G|>0\f$): the softened Coulomb / core-charge tail, the
    //! Fourier transform of the atom's long-range potential excluding \f$1/\Omega\f$ + the structure-factor
    //! phase (geometry, applied by the basis).  A pure Coulomb / Gaussian-smeared nucleus is ALL long, so this
    //! is the primary each model supplies. [energy x volume]
    virtual double FormFactorLong  (int Z, double G2) const=0;
    //! \brief SHORT-range (compact poly \f$\times\f$ Gaussian, no Coulomb tail) remainder.  Default: none (only
    //! a pseudopotential with a poly-Gaussian core, e.g. HGH, has one). [energy x volume]
    virtual double FormFactorShort (int Z, double G2) const {return 0.0;}
    //! \brief The full form factor \f$v(Z,|G|^2)=\f$ long + short.  Provided by the base -- a model supplies the
    //! pieces, never the sum.  UNIT-TEST-ONLY convenience: since the CP2K local-PP split (doc/GPWPlan.md 0e-PP)
    //! the production KS assembly consumes only the pieces (\c FormFactorLong -> the Hartree Poisson,
    //! \c FormFactorShort -> the external term); nothing in \c src/ calls the sum -- only the tests, as an oracle.
    virtual double FormFactor      (int Z, double G2) const {return FormFactorLong(Z,G2)+FormFactorShort(Z,G2);}

    //! \brief The \f$G\to0\f$ alignment of the LONG part: the finite \f$G\to0\f$ limit of \f$v_{long}(G)+4\pi
    //! Z/G^2\f$ (the softened-Coulomb remainder \f$\int[V_{long}+Z/r]\f$).  The "\f$\alpha\f$" a plane-wave
    //! total energy needs: the \f$G=0\f$ electron-ion energy is \f$(N_{el}/\Omega)\sum_a\alpha_a\f$, the
    //! uniform shift dropped with the \f$G=0\f$ potential.  Default 0 (a pure Coulomb tail has no remainder).
    virtual double FormFactorG0Long (int Z) const {return 0.0;}
    //! \brief The \f$G\to0\f$ alignment of the SHORT part.  Default: none.
    virtual double FormFactorG0Short(int Z) const {return 0.0;}
    //! \brief The full G=0 alignment \f$\alpha=\f$ long + short.  Provided by the base (model supplies pieces).
    //! UNIT-TEST-ONLY convenience (as \c FormFactor): production splits the alignment too -- \c FormFactorG0Long
    //! travels with \c FormFactorLong into the Hartree term, \c FormFactorG0Short stays with the external term.
    virtual double FormFactorG0     (int Z) const {return FormFactorG0Long(Z)+FormFactorG0Short(Z);}
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

//! One term \f$c\,r^{2n}\,e^{-\alpha r^2}\f$ of the SHORT-range local potential's closed Gaussian expansion
//! (\f$l=0\f$: a pure radial \f$\times\f$ even polynomial).  Structurally the \c SeparablePotential sibling
//! \c RadialGaussian, kept a distinct type so \c LocalPotential names no \c SeparablePotential module edge.
struct LocalGaussianTerm { double c; int n; double alpha; };

//! \brief CAPABILITY face: the SHORT-range part of a local pseudopotential in CLOSED GAUSSIAN form,
//! \f$V_{short}(Z,r)=\sum_t c_t\,r^{2n_t}\,e^{-\alpha_t r^2}\f$ -- i.e. the poly \f$\times\f$ Gaussian core
//! expressed exactly (not fitted) as radial-Gaussian terms.  A consumer holding a Gaussian orbital basis can
//! then assemble \f$\langle\chi_i|V_{short}|\chi_j\rangle\f$ ANALYTICALLY (a 3-centre Gaussian overlap,
//! \c Overlap3C) instead of a grid sweep -- the local sibling of \c SeparablePotential_Gaussian::BetaGaussian.
//! The LONG part (the erf-softened Coulomb tail) is NOT here: it folds into the G-space Poisson (doc/GPWPlan.md
//! 0e-PP).  Optional: reached by abstract->abstract cross-cast from \c LocalPotential; a model whose short part
//! is not Gaussian (or absent) simply does not implement it (the consumer keeps its grid path).
class LocalPotential_Gaussian : public virtual LocalPotential_Q
{
public:
    //! The closed Gaussian expansion of the SHORT-range part: \f$V_{short}(Z,r)=\sum_t c_t r^{2n_t}e^{-\alpha_t r^2}\f$.
    virtual std::vector<LocalGaussianTerm> ShortRangeGaussian(int Z) const=0;
};

//! \brief Bare nuclear Coulomb \f$v(G) = -4\pi Z/G^2\f$.  Physically exact but the 1s cusp makes the
//! plane-wave energy converge very slowly in \f$E_{cut}\f$.
class BareCoulomb : public LocalPotential
{
public:
    virtual double FormFactorLong(int Z, double G2) const {return -FourPi*Z/G2;}   // pure Coulomb: all long
    virtual double Vloc          (int Z, double r)  const {return -double(Z)/r;}    // -Z/r (singular at r=0)
};

//! \brief Gaussian-smeared nucleus (rung-1 "local pseudopotential"): the point charge is spread into a
//! Gaussian of width \f$\sigma\f$, giving \f$v(G) = -4\pi Z\,e^{-\sigma^2 G^2/2}/G^2\f$.  The smooth
//! charge removes the cusp, so the plane-wave energy converges rapidly with \f$E_{cut}\f$;
//! \f$\sigma\to 0\f$ recovers BareCoulomb.
class GaussianSmearedNucleus : public LocalPotential
{
public:
    explicit GaussianSmearedNucleus(double sigma) : itsSigma(sigma) {}
    //! A Gaussian core charge is a pure long-range (softened-Coulomb) field, no short remainder.
    virtual double FormFactorLong(int Z, double G2) const
    {
        return -FourPi*Z*std::exp(-0.5*itsSigma*itsSigma*G2)/G2;
    }
    //! \f$v(G)+4\pi Z/G^2 = 4\pi Z(1-e^{-\sigma^2G^2/2})/G^2 \to 2\pi Z\sigma^2\f$ as \f$G\to0\f$.
    virtual double FormFactorG0Long(int Z) const {return 2*Pi*Z*itsSigma*itsSigma;}
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
class HGH_LocalPotential : public LocalPotential, public virtual LocalPotential_Gaussian
{
public:
    //! \a Zion = valence charge, \a rloc = local radius, \a c = {C1,C2,C3,C4} polynomial coefficients.
    //! Real per-element parameters come from the GTH database via GetGTH (GTH_Potentials.C), not
    //! hardcoded factories.
    HGH_LocalPotential(double Zion, double rloc, double c1, double c2, double c3=0.0, double c4=0.0)
        : itsZion(Zion), itsRloc(rloc), itsC1(c1), itsC2(c2), itsC3(c3), itsC4(c4) {}

    //! \f$V_{long}(G)=-4\pi Z_{ion}\,e^{-G^2r_{loc}^2/2}/G^2\f$ -- the softened \f$-Z_{ion}/r\f$ tail (a
    //! Gaussian core charge of width \f$r_{loc}\f$).  Folded into the Hartree Poisson (doc/GPWPlan.md 0e-PP).
    virtual double FormFactorLong(int /*Z*/, double G2) const override
    {
        double t=G2*itsRloc*itsRloc;                          // (G r_loc)^2
        return -FourPi*itsZion/G2 * std::exp(-0.5*t);         // softened -Z_ion/r tail
    }
    //! \f$V_{short}(G)=(2\pi)^{3/2}r_{loc}^3\,e^{-G^2r_{loc}^2/2}\,\mathrm{poly}(t)\f$ -- the compact
    //! poly \f$\times\f$ Gaussian core (no Coulomb tail); the external term's remainder.
    virtual double FormFactorShort(int /*Z*/, double G2) const override
    {
        double t=G2*itsRloc*itsRloc;                          // (G r_loc)^2
        double g=std::exp(-0.5*t);
        double poly=itsC1 + itsC2*(3-t) + itsC3*(15-10*t+t*t) + itsC4*(105-105*t+21*t*t-t*t*t);
        double twopi32=std::pow(2*Pi, 1.5);
        return twopi32*itsRloc*itsRloc*itsRloc * g * poly;
    }
    // FormFactor == FormFactorLong + FormFactorShort is provided by the base (LocalPotential_Q).
    //! Real space (the closed-form inverse FT documented above):
    //! \f$V_{loc}(r)=-\tfrac{Z_{ion}}{r}\mathrm{erf}(\tfrac{r}{\sqrt2 r_{loc}})
    //! + e^{-r^2/2r_{loc}^2}\sum_i C_i (r/r_{loc})^{2i-2}\f$.  Finite at \f$r=0\f$ (no nuclear cusp).
    virtual double Vloc(int /*Z*/, double r) const override
    {
        double x=r/itsRloc, x2=x*x;
        double gpoly=std::exp(-0.5*x2)*(itsC1 + itsC2*x2 + itsC3*x2*x2 + itsC4*x2*x2*x2);
        double a=std::sqrt(2.0)*itsRloc;                       // erf argument scale
        double erfterm = (r>1e-12) ? -itsZion*std::erf(r/a)/r
                                   : -itsZion*2.0/(std::sqrt(Pi)*a);   // erf(r/a)/r -> 2/(sqrt(pi) a)
        return erfterm + gpoly;
    }
    // The full G=0 alignment \f$\alpha=\int[V_{loc}+Z_{ion}/r]\,d^3r = 2\pi Z_{ion}r_{loc}^2
    // + (2\pi)^{3/2}r_{loc}^3(C_1+3C_2+15C_3+105C_4)\f$ is FormFactorG0Long + FormFactorG0Short (the softened
    // Coulomb leaves \f$2\pi Z r_{loc}^2\f$, the Gaussian-polynomial the moment sum) -- summed by the base.
    //! \f$G\to0\f$ alignment of \c FormFactorLong: \f$2\pi Z_{ion}r_{loc}^2\f$ (the softened Coulomb remainder).
    virtual double FormFactorG0Long(int /*Z*/) const override
    {
        return 2*Pi*itsZion*itsRloc*itsRloc;
    }
    //! \f$G\to0\f$ alignment of \c FormFactorShort: \f$(2\pi)^{3/2}r_{loc}^3(C_1+3C_2+15C_3+105C_4)\f$.
    virtual double FormFactorG0Short(int /*Z*/) const override
    {
        double twopi32=std::pow(2*Pi, 1.5);
        double moments=itsC1 + 3*itsC2 + 15*itsC3 + 105*itsC4;
        return twopi32*itsRloc*itsRloc*itsRloc*moments;
    }
    virtual double Zion(int /*Z*/) const override {return itsZion;}   // valence charge (Z-independent)

    //! The SHORT part \f$V_{short}(r)=e^{-r^2/2r_{loc}^2}\sum_{i=1}^4 C_i(r/r_{loc})^{2(i-1)}\f$ as closed
    //! Gaussian terms: term \f$j\!=\!i\!-\!1\f$ is \f$c_j r^{2j}e^{-\alpha r^2}\f$ with \f$c_j=C_{j+1}/r_{loc}^{2j}\f$,
    //! \f$\alpha=1/2r_{loc}^2\f$ (term-by-term identical to \c Vloc's Gaussian-polynomial part, exact -- not a fit).
    virtual std::vector<LocalGaussianTerm> ShortRangeGaussian(int /*Z*/) const override
    {
        const double alpha=0.5/(itsRloc*itsRloc);
        const double C[4]={itsC1,itsC2,itsC3,itsC4};
        std::vector<LocalGaussianTerm> terms;
        for (int j=0;j<4;j++)
            if (C[j]!=0.0) terms.push_back({ C[j]/std::pow(itsRloc, 2*j), j, alpha });
        return terms;
    }
private:
    double itsZion, itsRloc, itsC1, itsC2, itsC3, itsC4;
};

//! \brief A multi-species local potential: a router keyed by atomic number \a Z that forwards each
//! query to the per-species sub-model.  This is ALL that multi-species (ionic) crystals need on the
//! local side -- every LocalPotential method already takes \a Z (single-species HGH ignores it; this
//! USES it to dispatch), so the basis assembly (which calls FormFactor(a->itsZ,...) per atom) is
//! unchanged.  Hand one of these to the external term and NaF / CsI just work.
class MultiSpecies_LocalPotential : public LocalPotential, public virtual LocalPotential_Gaussian
{
public:
    //! Register species \a Z's local model (atomic number Z, e.g. 11 for Na -- the atoms' itsZ).
    void Add(int Z, std::shared_ptr<const LocalPotential> model) {itsByZ[Z]=std::move(model);}
    // Forward the PRIMARY pieces per species; the base sums them into FormFactor / FormFactorG0.
    virtual double FormFactorLong   (int Z, double G2) const override {return Get(Z).FormFactorLong(Z,G2);}
    virtual double FormFactorShort  (int Z, double G2) const override {return Get(Z).FormFactorShort(Z,G2);}
    virtual double FormFactorG0Long (int Z)            const override {return Get(Z).FormFactorG0Long(Z);}
    virtual double FormFactorG0Short(int Z)            const override {return Get(Z).FormFactorG0Short(Z);}
    virtual double Vloc             (int Z, double r)  const override {return Get(Z).Vloc(Z,r);}
    virtual double Zion             (int Z)            const override {return Get(Z).Zion(Z);}
    //! The closed-Gaussian short view: cross-cast the sub-model to its Gaussian face (sanctioned
    //! abstract->abstract) -- mirrors \c MultiSpecies_SeparablePotential::BetaGaussian.
    virtual std::vector<LocalGaussianTerm> ShortRangeGaussian(int Z) const override
    {
        const auto* gface=dynamic_cast<const LocalPotential_Gaussian*>(&Get(Z));
        assert(gface && "MultiSpecies_LocalPotential::ShortRangeGaussian: sub-model has no closed-Gaussian view");
        return gface->ShortRangeGaussian(Z);
    }
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
