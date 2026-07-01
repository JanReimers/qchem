// File: Hamiltonian/Internal/Libxc_LDA.C  One LDA functional (exchange OR correlation) from libxc.
module;
#include <iosfwd>
#include <src/xc.h>
export module qchem.Hamiltonian.Internal.Libxc_LDA;
export import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Types;     // rvec3_t
import qchem.Vector3D;

export namespace qchem::Hamiltonian
{

//! A SINGLE LDA functional from libxc, selected by id (exchange OR correlation; see libxc.gitlab.io).
//! Provides BOTH the potential v (GetVxc -> xc_lda_vxc) and the energy density per particle eps_xc
//! (GetEpsXc -> xc_lda_exc).  Exposing eps_xc is the point: a CORRELATION instance feeds FittedVcorr the
//! correct E_c = integral eps_c rho rather than the 3/4 exchange-virial shortcut (which is exact only for
//! exchange).  Pair Libxc_LDA(1) [Dirac exchange] via FittedVxc with Libxc_LDA(corrId) via FittedVcorr to
//! assemble a full libxc LDA Hamiltonian -- the libxc analogue of SlaterExchange + VWN_Correlation.
//!
//! UNPOLARIZED-ONLY by construction: the scalar single-density face (GetVxc/GetEpsXc, one rho) cannot
//! supply the two channels an XC_POLARIZED libxc functional needs, so the functional always inits
//! XC_UNPOLARIZED and the single-density calls are always the correct libxc contract.  Polarized LDA is the
//! spin-native VWN5 path (Ham_DFTcorr_P / FittedVcorrPol), reached via XC::DiracVWN.
class Libxc_LDA : public ExFunctional
{
public:
    Libxc_LDA(int id);
    ~Libxc_LDA();   // xc_func_end (the matching teardown was missing before -> a libxc leak)

    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;
    virtual double  GetVxc  (double rho) const;   //!< v_xc(rho)         via xc_lda_vxc
    virtual double  GetEpsXc(double rho) const;   //!< eps_xc(rho)/part  via xc_lda_exc (overrides the 3/4 default)

    virtual std::ostream& Write(std::ostream&) const;

private:
    xc_func_type itsFunc;
};

} //namespace
