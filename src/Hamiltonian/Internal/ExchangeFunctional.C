// File: ExchangeFunctional.C   Exchange potential for DFT.
export module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Streamable;
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
};

} //namespace