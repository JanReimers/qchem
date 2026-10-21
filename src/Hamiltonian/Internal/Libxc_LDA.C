// File: Hamiltonian/Internal/Libxc_LDA.C  One LDA functional (exchange OR correlation) from libxc.
module;
#include <iosfwd>
#include <src/xc.h>
export module qchem.Hamiltonian.Internal.Libxc_LDA;
export import qchem.Hamiltonian.Internal.ExFunctional;

export namespace qchem::Hamiltonian
{

//! A SINGLE LDA functional from libxc, selected by id (exchange OR correlation; see libxc.gitlab.io).
//! Provides BOTH the potential v (GetVxc -> xc_lda_vxc) and the energy density per particle eps_xc
//! (GetEpsXc -> xc_lda_exc).  Exposing eps_xc is the point: a CORRELATION instance feeds FittedVcorr the
//! correct E_c = integral eps_c rho rather than the 3/4 exchange-virial shortcut (which is exact only for
//! exchange).  Pair Libxc_LDA(1) [Dirac exchange] via FittedVxc with Libxc_LDA(corrId) via FittedVcorr to
//! assemble a full libxc LDA Hamiltonian -- the libxc analogue of SlaterExchange + VWN_Correlation.
//!
//! SPIN-NATIVE (2026-10-21; was unpolarized-only): the functional is initialised TWICE, XC_UNPOLARIZED for the
//! scalar face and XC_POLARIZED for the two-channel face, so the libxc contract each face calls is the one its
//! arguments fit.  Exact \f$\zeta=0\f$ (\f$\rho_\uparrow=\rho_\downarrow\f$, every unpolarized run) stays on the
//! scalar path, byte-identical to before -- the same rule VWN_Correlation follows.
class Libxc_LDA : public ExFunctional
{
public:
    Libxc_LDA(int id);
    ~Libxc_LDA();   // xc_func_end (the matching teardown was missing before -> a libxc leak)

    virtual double  GetVxc  (double rho) const;   //!< v_xc(rho)         via xc_lda_vxc
    virtual double  GetEpsXc(double rho) const;   //!< eps_xc(rho)/part  via xc_lda_exc (overrides the 3/4 default)
    //! The spin-native face.  An EXCHANGE id (libxc kind XC_EXCHANGE) is channel-separable: the exact
    //! spin-scaling relation \f$v^\sigma=v_x(2\rho_\sigma)\f$, \f$\varepsilon^\sigma=\varepsilon_x(2\rho_\sigma)\f$ -- the
    //! base-class face over the scalar libxc call, the SAME per-channel split SlaterExchange uses.  The split
    //! matters: the energy is \f$\sum_\sigma\mathrm{Tr}(D_\sigma\,\varepsilon^\sigma_{fit})\f$, so two splits with the same
    //! \f$\sum_\sigma\rho_\sigma\varepsilon^\sigma\f$ still differ at the fit level (0.15 mHa on polarized N with libxc's
    //! one-\f$\varepsilon\f$ split).  Any other id goes through libxc's XC_POLARIZED contract (\c xc_lda_exc_vxc on
    //! \f$[\rho_\uparrow,\rho_\downarrow]\f$) with libxc's per-particle \f$\varepsilon_{xc}\f$ for BOTH channels, as
    //! VWN_Correlation does.  A negative channel
    //! (band-limited \f$\tilde\rho\f$ ringing) is clamped to 0 on its own (libxc's rule, and SlaterExchange's).
    virtual double  GetVxc  (double up, double dn, const Spin& s) const override;
    virtual double  GetEpsXc(double up, double dn, const Spin& s) const override;

    virtual std::ostream& Write(std::ostream&) const;

private:
    //! libxc at one (up, dn) point: eps per particle and (v_up, v_dn).  Clamps a negative channel to 0.
    void EvalPolarized(double up, double dn, double& eps, double v[2]) const;

    xc_func_type itsFunc;   //!< XC_UNPOLARIZED: the scalar face and exact zeta=0
    xc_func_type itsPol;    //!< XC_POLARIZED:   the two-channel face
    bool itsSeparable;      //!< an exchange id: the two-channel face is the spin-scaling relation
};

} //namespace
