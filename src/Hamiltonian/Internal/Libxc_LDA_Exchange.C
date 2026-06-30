// File: Hamiltonian/Internal/Libxc_LDA_Exchange.C Any LDA exchange potential defined in libxc.
module;
#include <iosfwd>
#include <src/xc.h>
export module qchem.Hamiltonian.Internal.Libxc_LDA_Exchange;
export import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Types;     // rvec3_t (was reached transitively via the removed Spin import)
import qchem.Vector3D;

export namespace qchem::Hamiltonian
{

//! An LDA exchange+correlation potential from libxc: Dirac exchange (LDA_X, id 1) + the correlation
//! functional named by \a id (e.g. 7 = LDA_C_VWN), summed.  See libxc.gitlab.io/functionals.
//!
//! UNPOLARIZED-ONLY BY CONSTRUCTION: the evaluation face is the scalar GetVxc(double rho) -- a SINGLE
//! density.  A libxc functional initialized XC_POLARIZED expects two channels (rho_up,rho_down) in and
//! (v_up,v_down) out, which this scalar interface cannot supply.  Rather than accept a Spin it can't honor,
//! the class simply does not offer a polarized mode (compile-time honesty, [[feedback_compile_time_over_runtime]]):
//! it always inits XC_UNPOLARIZED, so the single-density xc_lda_vxc call is always correct.  Polarized LDA
//! is the spin-native VWN5 path (Ham_DFTcorr_P / FittedVcorrPol), reached via XC::DiracVWN.
class Libxc_LDA_Exchange
    : public  ExFunctional
{
public:
    Libxc_LDA_Exchange(int id, double Ne);   // See https://libxc.gitlab.io/functionals/libxc-7.0.0/

    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;
    virtual double  GetVxc(double ChargeDensity) const;


    virtual std::ostream& Write(std::ostream&) const;

private:
    double Ne; //# of electrons
    xc_func_type exchange,corr;
};

} //namespace
