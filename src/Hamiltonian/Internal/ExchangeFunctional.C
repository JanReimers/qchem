// File: ExchangeFunctional.C   Exchange potential for DFT.
module;
#include <cassert>
#include <memory>
#include <vector>
#include <ostream>
export module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Streamable;
export import qchem.Symmetry.Spin;   // Spin -- the channel argument of the spin-native face
export import qchem.Types;           // rvec_t -- the batch records

export namespace qchem::Hamiltonian
{

//! \brief THE BATCH INPUT of an XC evaluation, spin-native (GGA design note §2; user ruling 2026-11-20).  One entry per point, in the
//! caller's order.  \c up / \c dn are \f$\rho_\uparrow,\rho_\downarrow\f$; the \f$\sigma\f$ triple
//! \f$\sigma_{\sigma\tau}=\nabla\rho_\sigma\cdot\nabla\rho_\tau\f$ is filled only for a functional that \c ConsumesDensityGradient (empty
//! otherwise: an LDA pays nothing).  An unpolarized run is its \f$\rho_\uparrow=\rho_\downarrow=\rho/2\f$,
//! \f$\sigma_{\sigma\tau}=\sigma/4\f$ collapse, which the TERM supplies (pin 5): no functional carries a polarization flag.
struct XCInput
{
    rvec_t up, dn;                 //!< \f$\rho_\uparrow,\rho_\downarrow\f$
    rvec_t sigUU, sigUD, sigDD;    //!< \f$\sigma_{\uparrow\uparrow},\sigma_{\uparrow\downarrow},\sigma_{\downarrow\downarrow}\f$; empty for an LDA
    size_t size() const {return up.size();}
    bool   HasGradient() const {return sigUU.size()!=0;}
};

//! \brief THE BATCH OUTPUT: the energy density PER UNIT VOLUME \f$e\f$ (what a quadrature integrates, \f$E_{xc}=\sum_g w_g e_g\f$), its
//! derivatives \f$v_{\rho\sigma}=\partial e/\partial\rho_\sigma\f$ and, for a GGA, \f$v_{\sigma\sigma\tau}=\partial e/\partial\sigma_{\sigma\tau}\f$
//! (libxc's \c vsigma layout).  The \f$v_\sigma\f$ triple is empty when the input carried no gradient.  A composite SUMS outputs.
struct XCOutput
{
    rvec_t e;                      //!< energy density per unit volume, \f$\sum_\sigma\rho_\sigma\varepsilon^\sigma\f$
    rvec_t vUp, vDn;               //!< \f$\partial e/\partial\rho_\uparrow,\ \partial e/\partial\rho_\downarrow\f$
    rvec_t vsUU, vsUD, vsDD;       //!< \f$\partial e/\partial\sigma_{\uparrow\uparrow},\ \partial e/\partial\sigma_{\uparrow\downarrow},\ \partial e/\partial\sigma_{\downarrow\downarrow}\f$; empty for an LDA input
    //! Zeroed and sized for \a in (the gradient triple only when \a in carries one): what an \c Evaluate accumulates into.
    void Reset(const XCInput& in)
    {
        const size_t n=in.size(), m=in.HasGradient() ? n : 0;
        e=rvec_t(n,0.0); vUp=rvec_t(n,0.0); vDn=rvec_t(n,0.0);
        vsUU=rvec_t(m,0.0); vsUD=rvec_t(m,0.0); vsDD=rvec_t(m,0.0);
    }
};

//! \brief An XC functional as a pure VALUE face: \f$\rho\mapsto v_{xc}\f$ and \f$\rho\mapsto
//! \varepsilon_{xc}\f$.  Data-free, like every other abstract interface in the project.
//!
//! It used to ALSO be a \c ScalarFunction<double> -- a FIELD face \f$r\mapsto v_{xc}(\rho(r))\f$ -- which
//! forced it to carry a \c const rChargeDensity* and a \c bool isPolarized as protected data, and to have
//! that pointer injected AFTER construction by \c InsertChargeDensity (V1.13's hidden-init landmine: on the
//! PW/GPW path it was simply never set, so \c op()(r) would have dereferenced null).  Every implementation
//! of that face was literally \c GetVxc((*itsChargeDensity)(r)), i.e. the value face composed with a density
//! the object should never have owned.
//!
//! A field is now built where it is USED, by an adapter holding (functional, density) as CONSTRUCTOR
//! arguments -- \c VxcDensity / \c EpsXcDensity in Imp/FittedVxc.C, \c PolVcDensity / \c PolEpsCDensity in
//! Imp/FittedVxc.C, the sampled fields in Imp/PWTerms_XC.C.  So the density arrives as an argument, not as
//! latched state, and the "which density am I attached to?" question cannot be answered wrongly.
//! (\c isPolarized went with it, and so -- V1.37 step 3 -- did \c SlaterExchange's Spin tag: polarization
//! is not a property of a functional at all; the TERM supplies per-channel densities to the spin-native
//! face below, and the imposed subgroup decides what those densities are.)
class ExFunctional
    : public virtual Streamable
{
public:
    //! \name The \f$\zeta=0\f$ scalar face: the TOTAL density of a closed shell.
    //! The historical single-density formulation, and still the primitive a CHANNEL-SEPARABLE functional
    //! (exchange) is written in.  Both are pure: \f$\varepsilon_{xc}=\tfrac34 v_x\f$ is exact for Dirac/Slater exchange
    //! only, so a default would be a silently wrong answer for any other functional (OOD §I.1).
    //!@{
    virtual double GetVxc  (double rho) const=0;
    virtual double GetEpsXc(double rho) const=0;
    //!@}

    //! \name THE SPIN-NATIVE FACE -- what every term consumes (V1.37 step 3, 2026-09-14).
    //!
    //! \f$v_{xc}^\sigma(\rho_\uparrow,\rho_\downarrow)\f$ and the per-channel energy per particle
    //! \f$\varepsilon^\sigma(\rho_\uparrow,\rho_\downarrow)\f$, defined so that
    //! \f$E_{xc}=\sum_\sigma\int\rho_\sigma\,\varepsilon^\sigma\f$ -- the form a per-spin density
    //! CONTRACTS against (\f$\mathrm{Tr}(D_\sigma\langle i|\varepsilon^\sigma|j\rangle)\f$), and the form a
    //! composite of functionals SUMS.  Spin-native is THE formulation; an unpolarized run is its
    //! \f$\rho_\uparrow=\rho_\downarrow=\rho/2\f$ collapse, which the TERM supplies (the folded doublet
    //! block hands over \f$\rho/2\f$ per channel) -- no functional carries a polarization flag any more.
    //!
    //! The DEFAULTS are the CHANNEL-SEPARABLE (exchange) rule, \f$v^\sigma=v_x(\rho_\sigma)\f$ with
    //! \f$v_x^\sigma(\rho_\sigma)=v_x^{\zeta=0}(2\rho_\sigma)\f$ (spin scaling: the closed-shell formula at
    //! the density a channel would have if doubled).  At \f$\zeta=0\f$ they reproduce the scalar face
    //! BIT FOR BIT (\f$2\cdot\tfrac12\rho=\rho\f$ exactly).  A correlation functional -- which COUPLES the
    //! channels through \f$r_s\f$ and \f$\zeta\f$ -- overrides both; a scalar-only wrapper (libxc)
    //! overrides them to refuse \f$\zeta\neq0\f$.
    //!@{
    virtual double GetVxc  (double up, double dn, const Spin& s) const {return GetVxc  (2.0*(s==Spin::Down ? dn : up));}
    virtual double GetEpsXc(double up, double dn, const Spin& s) const {return GetEpsXc(2.0*(s==Spin::Down ? dn : up));}
    //! The XC energy density PER UNIT VOLUME, \f$e=\sum_\sigma\rho_\sigma\varepsilon^\sigma\f$ -- what a
    //! quadrature integrates.  Per volume, not per particle, because per-particle forms do not share a
    //! denominator across a composite (\f$\rho_{tot}\f$ is zero in vacuum).  Derived, never overridden.
    double GetExcDensity(double up, double dn) const
    {return up*GetEpsXc(up,dn,Spin::Up) + dn*GetEpsXc(up,dn,Spin::Down);}
    //!@}

    //! \brief THE XC KERNEL \f$f_{xc}^{\sigma\sigma'}=\partial v^\sigma/\partial\rho_{\sigma'}\f$ at
    //! \f$(\rho_\uparrow,\rho_\downarrow)\f$ -- the second derivative of \f$E_{xc}\f$ a LINEAR RESPONSE needs
    //! (doc/LinearResponsePlan.md row H3).  Spin-native from the first line, like the face above (pin 5).
    //!
    //! THE DEFAULT IS A CENTRAL DIFFERENCE OF THIS FUNCTIONAL'S OWN spin-native \c GetVxc, pointwise, with a
    //! step relative to the local density -- so every functional answers, correlation included, with no
    //! second hand-derived formula to disagree with the first (VWN5's spin-polarised second derivatives are
    //! pages of algebra).  Its error is \f$O(h^2)\f$ in a smooth local function, ~1e-8 relative at the step
    //! used; an ANALYTIC override is the upgrade where one is short (Slater exchange has one).  Vanishing
    //! density (\f$\rho\le\f$ 1e-12) answers 0, the same guard every functional's \c GetVxc applies.
    virtual double GetFxc(double up, double dn, const Spin& s, const Spin& t) const
    {
        const double rho=up+dn;
        if (!(rho>1e-12)) return 0.0;
        const double h=1e-3*(t==Spin::Down ? dn : up)+1e-10*rho;   // relative to the channel being varied
        auto v=[&](double k){return t==Spin::Down ? GetVxc(up,dn+k*h,s) : GetVxc(up+k*h,dn,s);};
        return (8.0*(v(1)-v(-1))-(v(2)-v(-2)))/(12.0*h);           // four-point, O(h^4)
    }

    //! \brief How much denser the \f$v_{xc}\f$-fit grid must be than the wavefunction bandwidth, as a
    //! multiplier on the fit-basis energy cutoff (the CP2K \c REL_CUTOFF idea).
    //!
    //! The gradient enhancement of a GGA adds bandwidth to \f$v_{xc}\f$, so its fit grid must out-resolve the
    //! density's own {G}.  *How much* denser is a property of the FUNCTIONAL TYPE, which only the Hamiltonian
    //! side knows.  The functional distils that appetite to a scalar here; the basis stays functional-agnostic
    //! (handing an \c ExFunctional& to a \c qcBasisSet method would be a library cycle -- pass the NUMBER).
    //! LDA is band-limited to the density's own bandwidth \f$\Rightarrow 1\f$; a GGA overrides with ~1.5--2.
    virtual double GridCutoffFactor() const {return 1.0;}

    //! \name THE BATCH FACE (GGA design note §2; user ruling 2026-11-20, "Yes, all seven").
    //! One evaluation over many points, records in and out: what a GGA needs (libxc evaluates a batch in one call) and what the XC
    //! terms will move onto.  \c Evaluate ADDS this functional's contribution to \a out (sized by \c XCOutput::Reset), so a composite
    //! is a loop.  The DEFAULT loops the spin-native scalar face above, so every LDA (Slater, VWN, the libxc LDAs) answers with the
    //! numbers it gives today, point for point (a COMPOSITE's e associates its sum differently from its point face: roundoff, measured
    //! 1 ulp); a GGA overrides it (its scalar faces cannot answer without \f$\sigma\f$).
    //!@{
    virtual void Evaluate(const XCInput& in, XCOutput& out) const
    {
        for (size_t g=0;g<in.size();++g)
        {
            const double u=in.up[g], d=in.dn[g];
            out.e  [g]+=GetExcDensity(u,d);
            out.vUp[g]+=GetVxc(u,d,Spin::Up);
            out.vDn[g]+=GetVxc(u,d,Spin::Down);
        }
    }
    //! Does this functional read \f$\nabla\rho\f$?  What the CLIENT does with the answer: a term asks its sampler for the \f$\sigma\f$
    //! triple (and gathers the \f$\mathbf w\f$ field) only when this is true, so an LDA pays nothing; a policy refuses a route that
    //! cannot supply the gradient.
    virtual bool ConsumesDensityGradient() const {return false;}
    //!@}
};

// (The separate SpinCorrelation capability face is GONE -- V1.37 step 3.  It existed because only
//  correlation had a two-channel formulation while exchange rode a Spin-TAGGED functional and a bool on the
//  term; with spin-native as THE face on every functional, a term needs one face and no cross-cast.)

//! \brief A SUM OF FUNCTIONALS behaving as ONE (user, 2026-09-04).
//!
//! ⛔ WHY IT EXISTS, AND IT IS NOT TIDINESS.  Exchange and correlation were separate Hamiltonian TERMS,
//! and each term does its own real-space GATHER of its potential onto the basis.  The gather is LINEAR in
//! the field, so
//! \f$\langle i|v_x|j\rangle+\langle i|v_c|j\rangle=\langle i|(v_x+v_c)|j\rangle\f$: summing the two
//! POTENTIALS pointwise and gathering ONCE is the same operator for half the work.  Measured on the MnO
//! parity row (2026-09-04): 16 of 23 integrate-back calls per 3 SCF iterations were this split, against a
//! gather that is 71% of the run (doc/OpenWork.md bin 1).
//!
//! ⚠ NOT BIT-IDENTICAL: \c gather(a)+gather(b) and \c gather(a+b) differ in the last bits, so pinned
//! energies move at roundoff scale.  The OPERATOR is unchanged.
//!
//! EVERY part answers the same spin-native face with its own rule (channel-separable exchange reads its
//! own channel; correlation couples both), so the sum is a plain loop and adding a third functional is a
//! \c push_back, not a new term and not a new branch.
class CompositeExFunctional
    : public virtual ExFunctional
{
public:
    //! \a parts must be non-empty; each is kept alive by the composite.
    explicit CompositeExFunctional(std::vector<std::shared_ptr<ExFunctional>> parts) : itsParts(std::move(parts))
    {
        assert(!itsParts.empty() && "CompositeExFunctional: a sum of no functionals is not a functional");
    }

    // ---- the scalar (zeta=0) face: a plain sum -------------------------------------------------------
    virtual double GetVxc(double rho) const
    {double v=0.0; for (const auto& p : itsParts) v+=p->GetVxc(rho); return v;}
    virtual double GetEpsXc(double rho) const
    {double e=0.0; for (const auto& p : itsParts) e+=p->GetEpsXc(rho); return e;}
    // ---- the spin-native face: a plain sum, each part answering with its own rule ---------------------
    virtual double GetVxc(double up, double dn, const Spin& s) const
    {double v=0.0; for (const auto& p : itsParts) v+=p->GetVxc(up,dn,s); return v;}
    virtual double GetEpsXc(double up, double dn, const Spin& s) const
    {double e=0.0; for (const auto& p : itsParts) e+=p->GetEpsXc(up,dn,s); return e;}
    //! Each part answers with its OWN kernel (an analytic override where it has one) -- never a difference of the sum.
    virtual double GetFxc(double up, double dn, const Spin& s, const Spin& t) const
    {double f=0.0; for (const auto& p : itsParts) f+=p->GetFxc(up,dn,s,t); return f;}
    //! The batch: each part adds its own contribution (an LDA part leaves the \f$v_\sigma\f$ triple untouched, i.e. zero).
    virtual void Evaluate(const XCInput& in, XCOutput& out) const
    {for (const auto& p : itsParts) p->Evaluate(in,out);}
    //! Any part that reads \f$\nabla\rho\f$ makes the sum read it.
    virtual bool ConsumesDensityGradient() const
    {for (const auto& p : itsParts) if (p->ConsumesDensityGradient()) return true; return false;}
    //! The DENSEST part wins: a GGA in the mix sets the fit grid for the whole sum.
    virtual double GridCutoffFactor() const
    {double f=1.0; for (const auto& p : itsParts) f=max(f,p->GridCutoffFactor()); return f;}

    std::ostream& Write(std::ostream& os) const
    {for (const auto& p : itsParts) p->Write(os); return os;}

private:
    std::vector<std::shared_ptr<ExFunctional>> itsParts;
};

} //namespace