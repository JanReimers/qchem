// File: Hamiltonian/Internal/Libxc_GGA.C  One GGA functional (exchange OR correlation) from libxc.
module;
#include <iosfwd>
#include <src/xc.h>
export module qchem.Hamiltonian.Internal.Libxc_GGA;
export import qchem.Hamiltonian.Internal.ExFunctional;

export namespace qchem::Hamiltonian
{

//! \brief A SINGLE GGA functional from libxc, selected by id (exchange OR correlation; GGA design note §2, user ruling 2026-11-20).
//!
//! It answers ONLY the batch face: \f$e(\rho_\uparrow,\rho_\downarrow,\sigma_{\uparrow\uparrow},\sigma_{\uparrow\downarrow},\sigma_{\downarrow\downarrow})\f$
//! and its first derivatives, one \c xc_gga_exc_vxc call per batch under libxc's XC_POLARIZED contract.  The unpolarized run is the
//! term's \f$\rho/2\f$, \f$\sigma/4\f$ collapse fed to the same call (pin 5): one code path for both.  The scalar and spin-native POINT
//! faces cannot answer without \f$\sigma\f$, so they THROW -- a route that reaches them is a route the policy should have refused.
//!
//! libxc's own derivatives are consumed as they are (user ruling), each checked against a finite difference of libxc's \f$e\f$ by the
//! unit test (trust but verify).  libxc's density and \f$\sigma\f$ thresholds are this functional's own guards, printed by \c Write.
class Libxc_GGA : public ExFunctional
{
public:
    //! THROWS for an id libxc lacks, one it cannot initialise polarized, or one that is not a pure GGA (an LDA, a hybrid, a meta-GGA).
    explicit Libxc_GGA(int id);
    ~Libxc_GGA();

    //! \name The point faces: a GGA has no \f$\rho\f$-only value.  THROW (an internal invariant: the policy refuses every route that
    //! cannot supply \f$\nabla\rho\f$, so no term should get here).
    //!@{
    virtual double GetVxc  (double rho) const override;
    virtual double GetEpsXc(double rho) const override;
    virtual double GetVxc  (double up, double dn, const Spin& s) const override;
    virtual double GetEpsXc(double up, double dn, const Spin& s) const override;
    virtual double GetFxc  (double up, double dn, const Spin& s, const Spin& t) const override;
    //!@}

    //! \f$e\f$, \f$v_\rho\f$, \f$v_\sigma\f$ for the whole batch.  THROWS if \a in carries no \f$\sigma\f$ triple.  A negative channel
    //! (band-limited ringing) is clamped to 0 on its own, the rule \c Libxc_LDA follows; \f$\sigma\f$ is passed as given (libxc
    //! guards it with its own threshold).
    virtual void Evaluate(const XCInput& in, XCOutput& out) const override;
    virtual bool ConsumesDensityGradient() const override {return true;}

    virtual std::ostream& Write(std::ostream&) const override;
    //! Is this id a WHOLE exchange-correlation functional (libxc kind XC_EXCHANGE_CORRELATION), rather than an exchange or a
    //! correlation part?
    bool IsExchangeCorrelation() const;

private:
    xc_func_type itsPol;   //!< XC_POLARIZED: every evaluation (the unpolarized run is the term's collapse)
};

} //namespace
