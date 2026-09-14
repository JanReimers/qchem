// File: ExchangeFunctional.C   Exchange potential for DFT.
module;
#include <cassert>
#include <memory>
#include <vector>
#include <ostream>
export module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Streamable;
import qchem.Math;            // max (the composite's grid-cutoff reduction)
export import qchem.Symmetry.Spin;   // Spin -- the channel argument of the spin-native face

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
    //! (exchange) is written in.  \f$\varepsilon_{xc}\f$ defaults to the EXCHANGE virial \f$\tfrac34 v_x\f$
    //! (exact for Dirac/Slater exchange); CORRELATION functionals MUST override (differs ~15%).
    //!@{
    virtual double GetVxc  (double rho) const=0;
    virtual double GetEpsXc(double rho) const {return 0.75*GetVxc(rho);}
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
    //! The DENSEST part wins: a GGA in the mix sets the fit grid for the whole sum.
    virtual double GridCutoffFactor() const
    {double f=1.0; for (const auto& p : itsParts) f=max(f,p->GridCutoffFactor()); return f;}

    std::ostream& Write(std::ostream& os) const
    {for (const auto& p : itsParts) p->Write(os); return os;}

private:
    std::vector<std::shared_ptr<ExFunctional>> itsParts;
};

} //namespace