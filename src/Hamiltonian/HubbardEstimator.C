// File: Hamiltonian/HubbardEstimator.C  The U-ESTIMATOR face of DFT+U (increment 3, 2026-09-21): what a driver
// above the SCF consumes to obtain (U-bar, J-bar) per Hubbard manifold FROM THE CONVERGED ORBITALS -- ACBN0,
// Agapito/Curtarolo/Buongiorno Nardelli 2015 (doc/OpenWork.md step 5 increment 3, doc/Pins.md pin 23).
//
// WHY A FACE HERE.  The estimator needs the +U term's Löwdin projectors and the basis's bare on-site
// two-electron integrals, both of which live behind .Internal. modules a facade may not import; and it needs
// the ORBITALS (coefficients + occupations, per k and spin), which live ABOVE qcHamiltonian.  So the
// Hamiltonian hands out an estimator (tHamiltonian::MakeHubbardUEstimator, null when it carries no +U
// term), the driver FEEDS it the orbitals block by block, and asks for the estimate -- the same shape as the
// forward/adjoint pair: the data crosses the library boundary, the mechanism does not.
module;
#include <iosfwd>
#include <vector>
export module qchem.Hamiltonian.HubbardEstimator;
import qchem.BasisSet.Orbital_DFT_IBS;   // the block an orbital set lives on
import qchem.Symmetry.Irrep;             // Spin
import qchem.Types;

export namespace qchem::Hamiltonian
{

//! \brief One manifold's ACBN0 estimate.  Hartree throughout (the facade converts to eV at the edge).
struct HubbardEstimate
{
    size_t site=0;          //!< the manifold's site (cell atom order)
    int    l=2;             //!< its shell
    double Ubar=0;          //!< \f$\bar U\f$, eq 12: the renormalised on-site Coulomb average
    double Jbar=0;          //!< \f$\bar J\f$, eq 13: the renormalised on-site exchange average
    double Ueff() const {return Ubar-Jbar;}   //!< Dudarev's \f$U_{\rm eff}=\bar U-\bar J\f$, what the +U term takes
    rvec_t Nup, Ndn;        //!< \f$N^\sigma_m=\bar P^\sigma_{mm}\f$, the renormalised populations per function
    double chargeUp=0, chargeDn=0;   //!< their traces (the renormalised manifold charge per channel)
    double chargeUpBare=0, chargeDnBare=0;   //!< the unrenormalised manifold charge per channel (the pair-count populations)
    double UbarBare=0, JbarBare=0;   //!< the same averages with EVERY orbital weight 1 (the unrenormalised
                                     //!< shell averages of the manifold's actual occupation matrix) -- the
                                     //!< screening the renormalisation supplied is the ratio
};

//! \brief FEED me the occupied orbitals of every irrep block, then ask.  Occupations \a f are the block's
//! PHYSICAL ones (2 per spatial orbital under Spin::None, the folded doublet -- the estimator halves them
//! into both channels); \a w is the block's BZ weight; \a C holds one coefficient column per entry of \a f,
//! in the block's own AO basis (\c TOrbital::GetCoeff).
class HubbardUEstimator
{
public:
    virtual ~HubbardUEstimator() = default;
    virtual void Accumulate(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block, const Spin& s, double w,
                            const mat_t<double>& C, const rvec_t& f) = 0;
    virtual void Accumulate(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block, const Spin& s, double w,
                            const mat_t<dcmplx>& C, const rvec_t& f) = 0;
    //! One estimate per manifold, in the term's manifold order.  Throws before any orbital was fed.
    virtual std::vector<HubbardEstimate> Evaluate() const = 0;
    //! The estimates as a one-line report per manifold (eV): what the run banner and the probe print.
    virtual std::ostream& Write(std::ostream&) const = 0;
    //! WRITE the estimates into the term this estimator was built for: manifold \a M takes \f$U_{\rm eff}\f$ of
    //! \a e[M] for the next Fock build.  The paper's outer loop is Evaluate → Apply → re-converge → repeat;
    //! the facade's \c ConvergeHubbardU drives it.  Also resets the accumulated orbitals, so the next feed
    //! starts clean.
    virtual void Apply(const std::vector<HubbardEstimate>& e) = 0;
};

} // namespace
