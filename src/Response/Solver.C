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
//! ROWS are every MEASURED channel I (all of the probe's), COLUMNS the PERTURBED channels J -- which are an input
//! (doc/LinearResponsePlan.md §3d Q10): an inverse-response U depends on which channels the inverted matrix
//! spans (hp.x NiO: 5.267 eV perturbing Ni 3d, 5.434 eV adding O 2p), so the set is a statement of WHICH U.
struct SelfConsistentResponse
{
    std::vector<std::string> labels;    //!< every measured channel (the rows)
    std::vector<size_t> perturbed;      //!< the perturbed channels J, indices into \c labels (the columns)
    cmat_t              chi0;           //!< \f$\langle A^I,\mathcal R_0A^J\rangle\f$ -- no kernel.  I x J
    cmat_t              chi;            //!< \f$\langle A^I,\delta D^J\rangle\f$, \f$(1-\mathcal R_0\mathcal K)\delta D^J=\mathcal R_0A^J\f$.  I x J
    double              gap=0;          //!< the response gap the reference was gated on (E1)
    double              noise=0;        //!< the eigenvalue noise it was gated against (NaN = unmeasured)
    std::vector<double> residual;       //!< per PERTURBED channel: the Krylov relative residual reached
    std::vector<size_t> iterations;     //!< per PERTURBED channel: kernel applications
    //! The square J x J blocks (rows restricted to the perturbed set): what an inverse-response U inverts.
    cmat_t Chi0JJ() const;
    cmat_t ChiJJ () const;
    std::ostream& Write(std::ostream&) const;
};

//! \brief \f$\mathcal K\f$ applied to an orbital-basis δD and brought back to the orbital basis:
//! ToAO -> InducedFock -> ToMO, with δD RESCALED to unit max-norm on the way in and the result scaled back.
//! ⚠ WHY THE RESCALE (measured 2026-09-29, `ResponseKernel.GPW_Si_k211_KernelIsLinear_AtEveryScale`): the GPW
//! D-aware collocation screen has an ABSOLUTE tolerance, so the raw kernel is not homogeneous -- K[sδ]/s drifts
//! from K[δ] by 7e-5 at s = 1e-5 and 4% at s = 1e-8 (with the screen off: 1e-16).  A Krylov solve probes K with
//! unit-norm vectors of small components, and on NiO k222 its recurrence estimate ran 100x below the TRUE
//! residual.  K is linear, so the rescale changes nothing mathematically and makes the operator exactly
//! homogeneous; the screen then acts at the scale it was built for.
//! ▶ SINCE R3 STEP 3 THE PERIODIC KERNEL IS LINEAR OUTRIGHT at every q (the B1/B2 transition collocations take the
//! geometry-only screen; the same gate now asserts K[sδ]/s == K[δ] to 1e-12), so the rescale is a GUARD, not a fix:
//! it costs two scalings and keeps the operator homogeneous against any future term with an absolute tolerance.
template <class T> BlockPairs InducedFockMO(const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& K, const BlockPairs& dD,
                                            std::shared_ptr<const Symmetry::SelectionRule> rule);
//! \brief Solve the self-consistent response to each \a perturbed channel of \a probe under \a rule, measuring
//! every channel.  \a perturbed EMPTY = all of them (the molecular dipole, R2's single-manifold runs).  THROWS on
//! an index out of range.  FAILS on an inverted or unresolved coupled pair (E1) or a Krylov solve that does not
//! reach \a kp.tol (the residual is in the reason).
template <class T> Outcome<SelfConsistentResponse,ResponseFailure>
LinearResponse(const Reference& ref, const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& kernel,
               const ChannelProbe& probe, std::shared_ptr<const Symmetry::SelectionRule> rule, const KrylovParams& kp,
               std::vector<size_t> perturbed = {});

} // namespace
