// File: Response/Solver.C  The SELF-CONSISTENT linear response: (1 - R0 K) δD = R0 V by a Krylov solve
// (doc/LinearResponsePlan.md §1-§2, stage R1).
//
// THE FOUR OBJECTS MEET HERE, AND ONLY HERE (§2):
//   Reference  R0 : δF -> δD in the orbital basis (sum over states)       -- how independent particles respond
//   Kernel     K  : δD -> δF in the AO basis   (tHamiltonian's H1 face)   -- how the potential follows the density
//   Frame         : the AO <-> MO bridge between the two
//   Probe         : what we perturb (V_J) and what we read (the same operator, Adjoint/Forward)
// One Krylov step: unpack -> ToAO -> InducedFock -> ToMO -> ApplyR0 -> x - that.  No electronic-structure
// vocabulary leaks into the solver (M1), and no numerics leak into the Hamiltonian.
//
// A NON-CONVERGED SOLVE IS AN OUTCOME, NEVER A NUMBER (§7 trap 3): every channel's residual is reported, and a
// solve that stops above tolerance fails the whole response with the residual it reached.
module;
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>
export module qchem.Response.Solver;
export import qchem.Response.Probe;
export import qchem.Response.OrbitalFrame;
export import qchem.LASolver.Krylov;
export import qchem.Hamiltonian;   // ResponseKernel

export namespace qchem::Response
{

//! \brief The channel response at ONE selection rule (q = 0 / Invariant in R1): the bare \f$\chi_0\f$ and the
//! self-consistent \f$\chi\f$, with what each was gated and converged to.
struct SelfConsistentResponse
{
    std::vector<std::string> labels;
    cmat_t              chi0;           //!< \f$\langle A^I,\mathcal R_0A^J\rangle\f$ -- no kernel
    cmat_t              chi;            //!< \f$\langle A^I,\delta D^J\rangle\f$, \f$(1-\mathcal R_0\mathcal K)\delta D^J=\mathcal R_0A^J\f$
    double              gap=0;          //!< the response gap the reference was gated on (E1)
    double              noise=0;        //!< the eigenvalue noise it was gated against (NaN = unmeasured)
    std::vector<double> residual;       //!< per channel: the Krylov relative residual reached
    std::vector<size_t> iterations;     //!< per channel: kernel applications
    std::ostream& Write(std::ostream&) const;
};

//! \brief Solve the self-consistent response of every \a probe channel under \a rule.  FAILS on an inverted or
//! unresolved coupled pair (E1) or a Krylov solve that does not reach \a kp.tol (the residual is in the reason).
template <class T> Outcome<SelfConsistentResponse,ResponseFailure>
LinearResponse(const Reference& ref, const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& kernel,
               const ChannelProbe& probe, std::shared_ptr<const Symmetry::SelectionRule> rule, const KrylovParams& kp);

} // namespace
