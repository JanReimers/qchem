// File: Pseudopotential/LocalPotential.C  One-body LOCAL external potentials (plane-wave / pseudopotential).
//
// A local external potential is a SPECIES RADIAL FIELD (qchem.BasisSet.Orbital_PP_IBS): per nuclear species Z
// a radial function v_Z(r) with both spectral views, v_Z(r) and its form factor v~_Z(q^2).  The basis folds in
// the cell volume 1/Omega, the structure factor Sum_a e^{-iG.tau_a} and the G=0 handling
// (Orbital_PP_IBS::MakeSpeciesFieldMatrix); a LocalPotential supplies ONLY the one-species radial shape.
//
// THE CP2K LOCAL-PP SPLIT (doc/GPWPlan.md 0e-PP) is the model's natural structure: v = v_long + v_short, a
// softened-Coulomb / Gaussian core-charge tail plus a compact poly x Gaussian remainder.  A model supplies the
// PIECES (FormFactorLong/Short, VlocLong/Short, the two G->0 alignments) and this base answers the neutral
// face's FieldRange switch ONCE -- so the models read as the physics they are, and the basis never learns the
// word "pseudopotential".  This is the open/closed extension point for "lineage A": the bare Coulomb nucleus,
// the Gaussian-smeared nucleus and the analytic HGH/GTH local part are all just implementations of the pieces.
//
// V1.2 (2026-09-13): the abstract faces this module used to define (LocalPotential_Q/_R/_Gaussian) are GONE --
// they were structurally neutral already and only their names and home said "pseudopotential".  They live in
// qcBasisSet as SpeciesRadialField / SpeciesRadialField_Gaussian; this library implements them and depends on
// qcBasisSet, which no longer depends on it.
module;
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <vector>
#include <cassert>

export module qchem.Pseudopotential.LocalPotential;
export import qchem.BasisSet.Orbital_PP_IBS;   // SpeciesRadialField(+_Gaussian), FieldRange, Math::Gaussian
import qchem.Math; // FourPi, Pi

export namespace qchem::Pseudopotential
{

using BasisSet::FieldRange;
using BasisSet::SpeciesRadialField;
using BasisSet::SpeciesRadialField_Gaussian;

//! \brief A one-body local external potential: a \c SpeciesRadialField with the CP2K long/short PIECE
//! structure and its ion charge.  A model supplies the pieces; the neutral face's range switch is answered
//! here, once.  The full field is long + short, provided by this base -- a model never supplies the sum.
//!
//! The LONG part is the softened-Coulomb / Gaussian core-charge tail (folded into the G-space Poisson by the
//! periodic KS assembly); the SHORT part is the compact poly \f$\times\f$ Gaussian remainder (the external
//! term).  A pure Coulomb / Gaussian-smeared nucleus is ALL long, so \c FormFactorLong / \c VlocLong are the
//! primaries every model supplies and the short pieces default to none.
class LocalPotential : public virtual SpeciesRadialField
{
public:
    // --- the RECIPROCAL pieces: v~(Z,|G|^2), |G|>0, excluding 1/Omega + the structure-factor phase (geometry,
    //     applied by the basis). [energy x volume] ---
    virtual double FormFactorLong  (int Z, double G2) const=0;
    virtual double FormFactorShort (int Z, double G2) const {return 0.0;}
    // --- the REAL-space pieces: V(Z,r) in a.u. ---
    virtual double VlocLong  (int Z, double r) const=0;
    virtual double VlocShort (int Z, double r) const {return 0.0;}
    // --- the G->0 alignments: the finite G->0 limit of v~_long + 4 pi Z/G^2 (the softened-Coulomb remainder
    //     \f$\int[V_{long}+Z/r]\f$ -- the "alpha" a plane-wave total energy needs, the uniform shift dropped
    //     with the G=0 potential), and the short part's moment sum.  Default 0 (a pure Coulomb tail has none). ---
    virtual double FormFactorG0Long (int Z) const {return 0.0;}
    virtual double FormFactorG0Short(int Z) const {return 0.0;}

    //! \brief The ION charge this potential's \f$-Z_{ion}/r\f$ tail carries -- what the ion-ion (Ewald) sum
    //! needs.  Default = the (true) nuclear charge \a Z (all-electron: BareCoulomb, Gaussian-smeared); a
    //! pseudopotential overrides with its VALENCE charge (HGH: \f$Z_{ion}\f$).  A callback, not a getter: the
    //! answer depends on the model, not just on Z.  Consumed term-side only -- it is not part of the field.
    virtual double Zion(int Z) const {return double(Z);}
    //! Adapt Zion to the plain CALLBACK the ion-ion (Ewald) term consumes: Ewald (NuclearRepulsion) lives BELOW
    //! the pseudopotential layer in qcStructure, so it must take a neutral std::function.
    std::function<double(int)> ZionFn() const {return [this](int Z){return Zion(Z);};}

    // --- the neutral face, answered from the pieces ---
    virtual double ValueR(int Z, double r, FieldRange rng) const override
    {
        switch (rng)
        {
        case FieldRange::Long:  return VlocLong(Z,r);
        case FieldRange::Short: return VlocShort(Z,r);
        default:                return VlocLong(Z,r)+VlocShort(Z,r);
        }
    }
    virtual double ValueQ(int Z, double G2, FieldRange rng) const override
    {
        switch (rng)
        {
        case FieldRange::Long:  return FormFactorLong(Z,G2);
        case FieldRange::Short: return FormFactorShort(Z,G2);
        default:                return FormFactorLong(Z,G2)+FormFactorShort(Z,G2);
        }
    }
    virtual double CellMeanQ(int Z, FieldRange rng) const override
    {
        switch (rng)
        {
        case FieldRange::Long:  return FormFactorG0Long(Z);
        case FieldRange::Short: return FormFactorG0Short(Z);
        default:                return FormFactorG0Long(Z)+FormFactorG0Short(Z);
        }
    }
    // Conveniences that read as the physics (the sums the base owns).
    double FormFactor  (int Z, double G2) const {return ValueQ(Z,G2,FieldRange::Full);}
    double Vloc        (int Z, double r ) const {return ValueR(Z,r ,FieldRange::Full);}
    double FormFactorG0(int Z)            const {return CellMeanQ(Z,FieldRange::Full);}
};

//! \brief Bare nuclear Coulomb \f$v(G) = -4\pi Z/G^2\f$.  Physically exact but the 1s cusp makes the
//! plane-wave energy converge very slowly in \f$E_{cut}\f$.
class BareCoulomb : public LocalPotential
{
public:
    virtual double FormFactorLong(int Z, double G2) const {return -FourPi*Z/G2;}   // pure Coulomb: all long
    virtual double VlocLong      (int Z, double r)  const {return -double(Z)/r;}    // -Z/r (singular at r=0)
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
    virtual double VlocLong(int Z, double r) const
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
class HGH_LocalPotential : public LocalPotential, public virtual SpeciesRadialField_Gaussian
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
    // FormFactor == FormFactorLong + FormFactorShort is provided by the base.
    //! Real space, the closed-form inverse FT documented above, in the same two pieces: the LONG part
    //! \f$-\tfrac{Z_{ion}}{r}\mathrm{erf}(\tfrac{r}{\sqrt2 r_{loc}})\f$ (finite \f$-Z_{ion}\sqrt{2/\pi}/r_{loc}\f$ at
    //! \f$r=0\f$ -- no nuclear cusp) and the SHORT part \f$e^{-r^2/2r_{loc}^2}\sum_i C_i (r/r_{loc})^{2i-2}\f$.
    virtual double VlocLong(int /*Z*/, double r) const override
    {
        double a=std::sqrt(2.0)*itsRloc;                       // erf argument scale
        return (r>1e-12) ? -itsZion*std::erf(r/a)/r
                         : -itsZion*2.0/(std::sqrt(Pi)*a);     // erf(r/a)/r -> 2/(sqrt(pi) a)
    }
    virtual double VlocShort(int /*Z*/, double r) const override
    {
        double x=r/itsRloc, x2=x*x;
        return std::exp(-0.5*x2)*(itsC1 + itsC2*x2 + itsC3*x2*x2 + itsC4*x2*x2*x2);
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
    //! \f$\alpha=1/2r_{loc}^2\f$ (term-by-term identical to \c VlocShort, exact -- not a fit).  The LONG part
    //! (the erf-softened Coulomb tail) is not a Gaussian sum: it folds into the G-space Poisson instead, so
    //! asking for it (or Full) answers EMPTY -- "no closed form for this range", the consumer keeps its
    //! transform route.
    virtual std::vector<Math::Gaussian> AsGaussians(int /*Z*/, FieldRange rng) const override
    {
        std::vector<Math::Gaussian> terms;
        if (rng!=FieldRange::Short) return terms;
        const double alpha=0.5/(itsRloc*itsRloc);
        const double C[4]={itsC1,itsC2,itsC3,itsC4};
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
class MultiSpecies_LocalPotential : public LocalPotential, public virtual SpeciesRadialField_Gaussian
{
public:
    //! Register species \a Z's local model (atomic number Z, e.g. 11 for Na -- the atoms' itsZ).
    void Add(int Z, std::shared_ptr<const LocalPotential> model) {itsByZ[Z]=std::move(model);}
    // Forward the PIECES per species; the base sums them and answers the neutral face.
    virtual double FormFactorLong   (int Z, double G2) const override {return Get(Z).FormFactorLong(Z,G2);}
    virtual double FormFactorShort  (int Z, double G2) const override {return Get(Z).FormFactorShort(Z,G2);}
    virtual double FormFactorG0Long (int Z)            const override {return Get(Z).FormFactorG0Long(Z);}
    virtual double FormFactorG0Short(int Z)            const override {return Get(Z).FormFactorG0Short(Z);}
    virtual double VlocLong         (int Z, double r)  const override {return Get(Z).VlocLong(Z,r);}
    virtual double VlocShort        (int Z, double r)  const override {return Get(Z).VlocShort(Z,r);}
    virtual double Zion             (int Z)            const override {return Get(Z).Zion(Z);}
    //! The closed-Gaussian view: cross-cast the sub-model to its Gaussian face (sanctioned abstract->abstract)
    //! -- mirrors \c MultiSpecies_SeparablePotential::AsGaussians.
    virtual std::vector<Math::Gaussian> AsGaussians(int Z, FieldRange rng) const override
    {
        const auto* gface=dynamic_cast<const SpeciesRadialField_Gaussian*>(&Get(Z));
        assert(gface && "MultiSpecies_LocalPotential::AsGaussians: sub-model has no closed-Gaussian view");
        return gface->AsGaussians(Z,rng);
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
