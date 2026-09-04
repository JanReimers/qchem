// File: ExchangeFunctional.C   Exchange potential for DFT.
module;
#include <cassert>
#include <memory>
#include <vector>
#include <ostream>
export module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Streamable;
import qchem.Math;            // max (the composite's grid-cutoff reduction)
import qchem.Symmetry.Spin;   // Spin -- the requested channel of the spin-native correlation face

export namespace qchem::Hamiltonian
{

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
//! Imp/FittedVcorrPol.C, \c PWVxcField in Imp/PWTerms.C.  So the density arrives as an argument, not as
//! latched state, and the "which density am I attached to?" question cannot be answered wrongly.
//! (\c isPolarized went with it: it had NO setter caller, so it was permanently true, while \c GetVxc
//! actually branched on a DIFFERENT flag -- \c SlaterExchange::itsSpin.  Polarization is expressed by the
//! spin-native \c SpinCorrelation face below, and by that Spin, not by a bool on this face.)
class ExFunctional
    : public virtual Streamable
{
public:
    virtual double GetVxc(                double ChargeDensity) const=0;
    //! \brief Energy density per particle \f$\varepsilon_{xc}(\rho)\f$, so \f$E_{xc}=\int\varepsilon_{xc}\rho\,d^3r\f$.
    //! Default is the EXCHANGE virial \f$\varepsilon_x=\tfrac34 v_x\f$ (exact for Dirac/Slater exchange).
    //! CORRELATION functionals MUST override: \f$\varepsilon_c\neq\tfrac34 v_c\f$ (differs ~15%).
    virtual double GetEpsXc(              double ChargeDensity) const {return 0.75*GetVxc(ChargeDensity);}

    //! \brief How much denser the \f$v_{xc}\f$-fit grid must be than the wavefunction bandwidth, as a
    //! multiplier on the fit-basis energy cutoff (the CP2K \c REL_CUTOFF idea).
    //!
    //! The gradient enhancement of a GGA adds bandwidth to \f$v_{xc}\f$, so its fit grid must out-resolve the
    //! density's own {G}.  *How much* denser is a property of the FUNCTIONAL TYPE, which only the Hamiltonian
    //! side knows.  The functional distils that appetite to a scalar here; the basis stays functional-agnostic
    //! (handing an \c ExFunctional& to a \c qcBasisSet method would be a library cycle -- pass the NUMBER).
    //! LDA is band-limited to the density's own bandwidth \f$\Rightarrow 1\f$; a GGA overrides with ~1.5--2.
    virtual double GridCutoffFactor() const {return 1.0;}
};

//! \brief Spin-native correlation face (no data; an abstract capability mixin).
//!
//! Correlation does NOT separate by spin channel the way Slater exchange does: \f$v_c^\sigma\f$ and
//! \f$\varepsilon_c\f$ COUPLE both densities (through \f$r_s(\rho_\uparrow+\rho_\downarrow)\f$ and
//! \f$\zeta\f$), so they cannot be expressed through the single-density \c ExFunctional::GetVxc face that
//! channel-separable exchange uses.  A correlation functional that supports polarized DFT implements this
//! two-channel face; \c FittedVcorrPol consumes it.  Unpolarized is the \f$\rho_\uparrow=\rho_\downarrow\f$
//! collapse (so \c GetVc(h,h,s)==\c ExFunctional::GetVxc(2h)).
class SpinCorrelation
{
public:
    virtual ~SpinCorrelation() {}
    //! Correlation energy density per particle \f$\varepsilon_c(\rho_\uparrow,\rho_\downarrow)\f$.
    virtual double GetEpsC(double rhoUp, double rhoDown) const=0;
    //! Channel correlation potential \f$v_c^\sigma=\varepsilon_c+\rho\,\partial\varepsilon_c/\partial\rho_\sigma\f$.
    virtual double GetVc  (double rhoUp, double rhoDown, const Spin&) const=0;

    //! \brief The XC energy DENSITY **PER UNIT VOLUME**, \f$e(\rho_\uparrow,\rho_\downarrow)\f$, so that
    //! \f$E=\int e\,d^3r\f$.
    //!
    //! ⚠ PER VOLUME, NOT PER PARTICLE, and that is the whole reason it exists: a COMPOSITE of functionals
    //! sums energy densities, and the per-particle forms do NOT share a denominator -- exchange carries
    //! \f$\varepsilon_x(\rho_\sigma)\rho_\sigma\f$ summed over channels while correlation carries
    //! \f$\varepsilon_c\rho_{tot}\f$.  Adding them per particle would need a division by \f$\rho_{tot}\f$,
    //! which is ZERO in vacuum.  Per volume they just add.
    //! Default: the single-correlation-functional form, so every existing functional is unchanged.
    virtual double GetExcDensity(double rhoUp, double rhoDown) const
    {return GetEpsC(rhoUp,rhoDown)*(rhoUp+rhoDown);}
};

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
//! IT IMPLEMENTS BOTH FACES, because its parts may be either.  A part that is only an \c ExFunctional is
//! CHANNEL-SEPARABLE (exchange): through the two-channel face it contributes \f$v_x(\rho_\sigma)\f$,
//! reading its own channel and ignoring the other.  A part that is also a \c SpinCorrelation contributes
//! its coupled \f$v_c^\sigma(\rho_\uparrow,\rho_\downarrow)\f$.  ⇒ Adding a third functional is a
//! \c push_back, not a new term and not a new branch.
class CompositeExFunctional
    : public virtual ExFunctional
    , public virtual SpinCorrelation
{
public:
    //! \a parts must be non-empty; each is kept alive by the composite.
    //! \note THE CROSS-CAST IS RESOLVED ONCE, HERE.  \c GetVc / \c GetExcDensity run PER GRID POINT
    //! (~64k points x 2 channels per matrix build), and a \c dynamic_cast on that path would be a
    //! per-point RTTI lookup for a fact that is fixed at construction.  It is an abstract-to-abstract
    //! cast between two capability faces, which is the legitimate kind (CLAUDE.md) -- it just has no
    //! business in an inner loop.
    explicit CompositeExFunctional(std::vector<std::shared_ptr<ExFunctional>> parts)
    {
        assert(!parts.empty() && "CompositeExFunctional: a sum of no functionals is not a functional");
        for (auto& f : parts)
        {
            const SpinCorrelation* c=dynamic_cast<const SpinCorrelation*>(f.get());
            itsParts.push_back(Part{std::move(f), c});
        }
    }

    // ---- ExFunctional (the unpolarized face): a plain sum -------------------------------------------
    virtual double GetVxc(double rho) const
    {double v=0.0; for (const Part& p : itsParts) v+=p.f->GetVxc(rho); return v;}
    virtual double GetEpsXc(double rho) const
    {double e=0.0; for (const Part& p : itsParts) e+=p.f->GetEpsXc(rho); return e;}
    //! The DENSEST part wins: a GGA in the mix sets the fit grid for the whole sum.
    virtual double GridCutoffFactor() const
    {double f=1.0; for (const Part& p : itsParts) f=max(f,p.f->GridCutoffFactor()); return f;}

    // ---- SpinCorrelation (the spin-native face) ------------------------------------------------------
    virtual double GetVc(double up, double dn, const Spin& s) const
    {
        double v=0.0;
        for (const Part& p : itsParts)
            if (p.spin) v+=p.spin->GetVc(up,dn,s);
            else        v+=p.f->GetVxc(s==Spin::Down ? dn : up);   // channel-separable: its OWN channel
        return v;
    }
    //! \note Present because the face demands it; the TERM integrates \c GetExcDensity instead, which is
    //! the form that composes.  Dividing by \f$\rho_{tot}\f$ here would be undefined in vacuum, so this
    //! reports the per-particle value only where that is meaningful and 0 where it is not.
    virtual double GetEpsC(double up, double dn) const
    {const double r=up+dn; return r>0.0 ? GetExcDensity(up,dn)/r : 0.0;}
    virtual double GetExcDensity(double up, double dn) const
    {
        double e=0.0;
        for (const Part& p : itsParts)
            if (p.spin) e+=p.spin->GetExcDensity(up,dn);
            else        e+=p.f->GetEpsXc(up)*up + p.f->GetEpsXc(dn)*dn;   // E_x = Σ_σ ∫ ε_x(ρ_σ)ρ_σ
        return e;
    }

    std::ostream& Write(std::ostream& os) const
    {for (const Part& p : itsParts) p.f->Write(os); return os;}

private:
    //! One summand, with its two-channel capability resolved at construction.
    struct Part
    {
        std::shared_ptr<ExFunctional> f;            //!< owns the summand
        const SpinCorrelation*        spin=nullptr; //!< non-null iff \c f also has the coupled face
    };
    std::vector<Part> itsParts;
};

} //namespace