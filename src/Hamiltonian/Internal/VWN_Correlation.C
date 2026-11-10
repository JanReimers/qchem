// File: VWN_Correlation.C  Vosko-Wilk-Nusair (VWN5) LDA correlation functional.
module;
#include <ostream>
export module qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Symmetry.Spin;   // Spin (the requested correlation channel)

export namespace qchem::Hamiltonian
{

//! \brief VWN5 (Vosko-Wilk-Nusair, Can.J.Phys 58, 1200 (1980), functional V) LDA CORRELATION.
//!
//! SPIN-NATIVE: the primary formulation is the two-channel \f$\varepsilon_c(\rho_\uparrow,\rho_\downarrow)\f$
//! with the per-channel potentials \f$v_c^\sigma\f$ (the two-channel GetEpsXc / GetVxc).  The spin-unpolarized scalar face
//! (GetVxc / GetEpsXc) is the \f$\zeta=0\f$ collapse, kept as its own path because it is cheaper (one G, not
//! three); the two-channel faces route exact \f$\zeta=0\f$ through it.
//!
//! DERIVATIVES ARE ENZYME'S (2026-11-09, the OpenWork "Enzyme REQUIRED" row, family 2).  The functional is written
//! ONCE, as its energy per volume \f$e_c(\rho_\uparrow,\rho_\downarrow)=\rho\,\varepsilon_c\f$ (Imp/VWN_Correlation.C);
//! \f$v_c^\sigma=\partial e_c/\partial\rho_\sigma\f$ is its reverse-mode derivative and the kernel
//! \f$f_c^{\sigma\sigma'}=\partial v_c^\sigma/\partial\rho_{\sigma'}\f$ the forward-mode derivative of that.  The
//! hand \f$dG/dx\f$, \f$f'(\zeta)\f$ and \f$(r_s,\zeta)\f$ partials are gone, and so is the base class's
//! finite-difference kernel for VWN.  The FD gates stay (tests/LDA_XC.C) while EnzymeAD/Enzyme#3747 is open.
//!
//! Validated pointwise against libxc LDA_C_VWN (see LDA_XC_UT), both XC_UNPOLARIZED and XC_POLARIZED.
//! Provides BOTH the potential \f$v_c\f$ and the energy density \f$\varepsilon_c\f$ -- the latter is
//! essential because \f$E_c=\int\varepsilon_c\rho\f$ and \f$\varepsilon_c\neq\tfrac34 v_c\f$, so the
//! exchange-virial energy shortcut used for Slater exchange is invalid for correlation.
//!
//! Spin interpolation (VWN Eq. 4.4, the form used by libxc):
//! \f[ \varepsilon_c(r_s,\zeta)=\varepsilon_c^P + \alpha_c(r_s)\frac{f(\zeta)}{f''(0)}(1-\zeta^4)
//!                              + [\varepsilon_c^F-\varepsilon_c^P]\,f(\zeta)\,\zeta^4 \f]
//! with \f$f(\zeta)=\frac{(1+\zeta)^{4/3}+(1-\zeta)^{4/3}-2}{2^{4/3}-2}\f$ interpolating paramagnetic
//! (\f$\zeta=0\f$) to ferromagnetic (\f$\zeta=1\f$), and \f$\alpha_c\f$ the spin stiffness.
class VWN_Correlation : public ExFunctional
{
public:
    VWN_Correlation() {}

    // --- spin-unpolarized scalar face (the zeta=0 collapse) ---
    virtual double GetVxc  (double rho) const;
    virtual double GetEpsXc(double rho) const;

    // --- spin-native two-channel face (the primary formulation) ---
    //! Correlation energy density per particle \f$\varepsilon_c(\rho_\uparrow,\rho_\downarrow)\f$ -- the same
    //! for both channels (the Spin argument is the face's, not this functional's).  At EXACT \f$\zeta=0\f$
    //! it takes the scalar path.
    virtual double GetEpsXc(double rup, double rdn, const Spin&) const;
    //! Channel correlation potential \f$v_c^\sigma=\partial e_c/\partial\rho_\sigma\f$ (Enzyme).
    //! Note v_c^sigma COUPLES both channels (through r_s and zeta) -- it is NOT a function of rho_sigma alone.
    virtual double GetVxc(double rup, double rdn, const Spin& s) const;
    //! The kernel \f$f_c^{\sigma\sigma'}=\partial v_c^\sigma/\partial\rho_{\sigma'}\f$ (Enzyme, forward over reverse).
    //! At full polarization it diverges (\f$f''(\zeta)\propto(1\mp\zeta)^{-2/3}\f$), so the minority channel is
    //! floored at a stated \f$\zeta\f$ threshold first (Imp/VWN_Correlation.C, \c kZetaThreshold).
    virtual double GetFxc(double rup, double rdn, const Spin& s, const Spin& t) const override;

    virtual std::ostream& Write(std::ostream& os) const { return os << "VWN5"; }
};

} //namespace
