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

//! \brief The self-consistent channel response over a q-MESH (doc/LinearResponsePlan.md §3d, R3 step 4): one
//! \c SelfConsistentResponse per q (one ground state, one kernel; the q points loop INSIDE the response), and hp.x's
//! real-space quantities over the q-mesh's SUPERCELL -- the response to an ISOLATED perturbation, which is what a
//! supercell of that size at Γ gives (the supercell-equivalence gate).  At \f$N_q=1\f$ it is the q = 0 response of
//! this cell, and the supercell is the cell.
struct QMeshResponse
{
    ivec3_t                             Nq;
    std::vector<MeshShift>              q;      //!< the q-mesh, in \c Reference::QMesh order
    std::vector<SelfConsistentResponse> atQ;    //!< per q; every entry has the same labels and perturbed set
    //! The J x J (perturbed) blocks in real space, rows/columns \f$(J,\mathbf R)\f$ at \f$r\,n_J+c\f$ (\c ToRealSpace's layout).
    //! THROWS if the result is not real (time reversal): \f$\chi_0\f$ to 1e-8 relative (it is exact), \f$\chi\f$ to 1e4 x the
    //! worst Krylov residual -- a residual r bounds the solution error only up to the conditioning of
    //! \f$(1-\mathcal R_0\mathcal K)\f$, measured ~190 r on NiO k222 (2026-10-04: leak 1.6e-6 at r = 8.7e-9).  \a imagRatio:
    //! the leak measured (Write prints it).
    rmat_t Chi0RealSpace(double* imagRatio=nullptr) const;
    rmat_t ChiRealSpace (double* imagRatio=nullptr) const;
    //! \f$W=\chi_0^{-1}-\chi^{-1}\f$ over the WHOLE supercell's J x J real-space matrices (Hartree; \c ToRealSpace's layout).  Its
    //! home-cell diagonal is U, its home-cell off-diagonals the inter-site \f$V_{IJ}\f$ (Timrov et al. PRB 98, 085127 eq (2);
    //! PRB 103, 045141 for V) -- over the PERTURBED set only, so a V to a site outside it is not here.
    rmat_t Interaction() const;
    //! \f$U_c=W_{(c,\mathbf 0),(c,\mathbf 0)}\f$ (Hartree) per perturbed channel c: hp.x's definition.
    rvec_t U() const;
    //! \name THE ORBITAL-RESOLVED INVERSIONS (OpenWork §2 "Response: REPORT orbital-resolved channels", 2026-10-21).
    //! Inverting the full matrix over the resolved channels of ONE shell is the near-singular case (doc/Pins.md pin 24):
    //! Macke et al. JCTC 20, 4824 (2024) invert "pointwise".  Both of their readings, beside hp.x's full one:
    //!@{
    //! \f$W\f$ with \f$\chi_0\f$ and \f$\chi\f$'s off-diagonals between perturbed channels of one GROUP in one cell
    //! ZEROED before inverting (a group = one site's orbitals: inter-orbital off-diagonals dropped, inter-site kept;
    //! each diagonal \f$\chi_{cc}\f$ still carries the intra-shell screening).  \a group: per PERTURBED channel.
    rmat_t Interaction(const std::vector<size_t>& group) const;
    rvec_t U(const std::vector<size_t>& group) const;
    //! \f$U_c=1/\chi_{0,cc}-1/\chi_{cc}\f$ in the home cell -- a single channel perturbed and measured alone in the
    //! q-mesh's supercell: what a \f$\pm\alpha\f$ finite difference on one orbital (pw.x orbital_resolved) measures.
    rvec_t UScalar() const;
    //!@}
    //! Per MEASURED channel (the probe's): a caller's grouping key -- the facade's Hubbard MANIFOLD, so a resolved channel
    //! can be put back on its site.  EMPTY = every channel its own (a manifold-resolved response).
    std::vector<size_t> group;
    //! The supercell's real-space home-cell blocks, the time-reversal leak and U (eV, hp.x's unit).  (Each q's own table is
    //! printed by \c LinearResponse as that q finishes.)
    std::ostream& Write(std::ostream&) const;
};

//! \brief THE LINEAR-RESPONSE U ESTIMATOR (doc/LinearResponsePlan.md R4, H5): hp.x's \f$U=(\chi_0^{-1}-\chi^{-1})_{II}\f$ of a solved
//! q-mesh response, one estimate per PERTURBED channel, behind the same strategy face as ACBN0.  Its input is its own
//! construction business (D4): the facade solves the response (+U frozen) and hands it over with the term's channels, which name
//! the manifolds.  \a channels[J] must describe the probe's channel J (the probe is made from the same \c HubbardChannels).
class LR_HubbardU : public virtual Hamiltonian::HubbardUEstimator
{
public:
    LR_HubbardU(QMeshResponse r, std::vector<Hamiltonian::HubbardChannel> channels);
    virtual std::vector<Hamiltonian::HubbardEstimate> Evaluate() const override;
    virtual std::ostream& Write(std::ostream&) const override;
    const QMeshResponse& Response() const {return itsR;}
private:
    QMeshResponse                           itsR;
    std::vector<Hamiltonian::HubbardChannel> itsChannels;
};

//! \brief \c LinearResponse at every q of the \a Nq mesh (\c Reference::QMesh), each its own Krylov solve on the
//! SAME reference, frame and kernel.  FAILS on an incommensurate q-mesh, or with the first q whose response fails
//! (the q is named in the reason).
template <class T> Outcome<QMeshResponse,ResponseFailure>
LinearResponse(const Reference& ref, const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& kernel,
               const ChannelProbe& probe, ivec3_t Nq, const KrylovParams& kp, std::vector<size_t> perturbed = {});

} // namespace
