// File: Hamiltonian/HubbardEstimator.C  The U-ESTIMATOR faces of DFT+U: the STRATEGY face the outer U loop consumes
// (doc/LinearResponsePlan.md H5, ruling D4: ACBN0 or linear response, the facade picks), and ACBN0's own INPUT face --
// (U-bar, J-bar) per Hubbard manifold FROM THE CONVERGED ORBITALS, Agapito/Curtarolo/Buongiorno Nardelli 2015
// (doc/OpenWork.md step 5 increment 3, doc/Pins.md pin 23).
//
// H5 (2026-10-05): ESTIMATING a U and APPLYING it are two responsibilities.  An estimator answers "what U?"; the outer loop
// writes the answer through \c HubbardUTarget::SetU, whichever strategy produced it.  Each strategy's INPUT is its own
// construction business: ACBN0 is fed orbitals (\c ACBN0Estimator::Accumulate, by qcResponse's FeedOccupiedOrbitals); the
// linear-response estimator (qcResponse's LR_HubbardU) is built from a solved response.
//
// WHY A FACE HERE.  The estimator needs the +U term's Löwdin projectors and the basis's bare on-site
// two-electron integrals, both of which live behind .Internal. modules a facade may not import; and it needs
// the ORBITALS (coefficients + occupations, per k and spin), which live ABOVE qcHamiltonian.  So the
// Hamiltonian hands out ACBN0's estimator (tHamiltonian::MakeACBN0Estimator, null when it carries no +U
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

//! \brief One manifold's estimated U -- what ANY estimator hands the outer loop.  Hartree (the facade converts to eV at the edge).
struct HubbardEstimate
{
    size_t manifold=0;      //!< index into the +U term's manifolds (what \c HubbardUTarget::SetU takes)
    size_t site=0;          //!< the manifold's site (cell atom order) -- for the report
    int    l=2;             //!< its shell
    double U=0;             //!< the U the +U term should take (Dudarev's U_eff)
};

//! \brief THE STRATEGY FACE: what the outer U loop (\c SolidCalculation::ConvergeHubbardU) consumes.  One estimate per manifold
//! the strategy sets -- ACBN0 every manifold, linear response its perturbed (Hubbard-site) set.
class HubbardUEstimator
{
public:
    virtual ~HubbardUEstimator() = default;
    virtual std::vector<HubbardEstimate> Evaluate() const = 0;
    //! One line per manifold (eV): the strategy reports what it estimated, at its own activity (pin 17).
    virtual std::ostream& Write(std::ostream&) const = 0;
};

//! \brief ACBN0's DETAIL for one manifold (eqs 12-13) -- its own result, not the strategy face's (H5).
struct ACBN0Estimate
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

//! \brief ACBN0's INPUT face: FEED me the occupied orbitals of every irrep block, then ask.  Occupations \a f are the block's
//! PHYSICAL ones (2 per spatial orbital under Spin::None, the folded doublet -- the estimator halves them into both channels);
//! \a w is the block's BZ weight; \a C holds one coefficient column per entry of \a f, in the block's own AO basis
//! (\c TOrbital::GetCoeff).  \c Evaluate answers \f$U_{\rm eff}\f$ for every manifold; \c Details the whole ACBN0 picture.
class ACBN0Estimator : public virtual HubbardUEstimator
{
public:
    virtual void Accumulate(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block, const Spin& s, double w,
                            const mat_t<double>& C, const rvec_t& f) = 0;
    virtual void Accumulate(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block, const Spin& s, double w,
                            const mat_t<dcmplx>& C, const rvec_t& f) = 0;
    //! One per manifold, in the term's manifold order.  Throws before any orbital was fed.
    virtual std::vector<ACBN0Estimate> Details() const = 0;
};

} // namespace
