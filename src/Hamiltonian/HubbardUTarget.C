// File: Hamiltonian/HubbardUTarget.C  What a U-ESTIMATOR does TO the +U term, as a public capability
// (doc/LinearResponsePlan.md §3d, Q6 -- H5's HubbardUTarget, pulled forward for R3).
//
// WHY A FACE HERE.  The +U term lives behind the .Internal. Hubbard module, and the facade that runs a linear
// response sits above qcHamiltonian.  Linear-response U is DEFINED with V_Hub held at its ground-state value
// (Timrov et al. PRB 98, 085127, eq 20), so the facade must freeze +U for the length of a response solve --
// and the finite-difference LRT cross-check (Cococcioni & de Gironcoli 2005) must apply a static alpha*P to one
// manifold of a CONVERGED run and re-converge with +U frozen.  Both are writes TO the term; the projectors a
// response perturbs and measures through are the read-only sibling, HubbardChannels.  Two faces because they
// are two reasons to change (a probe never writes the term; an estimator never needs the amplitudes).
module;
#include <cstddef>
export module qchem.Hamiltonian.HubbardUTarget;

export namespace qchem::Hamiltonian
{
//! \brief The WRITES a U estimator makes to the run's +U term.  Handed out by \c tHamiltonian::GetHubbardUTarget
//! (null when the run has no +U term); keep the Hamiltonian alive while it is used.  Every write takes effect at
//! the NEXT Fock build, frozen or not.
class HubbardUTarget
{
public:
    virtual ~HubbardUTarget() = default;
    //! Set manifold \a M's \f$U\f$ (Hartree; a filled \c Uirrep is set to it throughout, shell-averaged) -- the
    //! outer loop's write (ACBN0 today, linear response at R4).
    virtual void SetU(size_t M, double U) = 0;
    //! \brief Set the STATIC perturbation \f$\alpha\hat P_M\f$ on manifold \a M (Hartree; QE's \c Hubbard_alpha):
    //! \f$W\mathrel{+}=\alpha\mathbb 1\f$, \f$E\mathrel{+}=\alpha\,{\rm Tr}\,n\f$.  The LR-cDFT perturbation: converge
    //! at \f$\pm\alpha\f$ and difference the occupation, and that IS \f$\chi=dn/d\alpha\f$.  0 removes it.
    virtual void SetPerturbation(size_t M, double alpha) = 0;
    //! \brief Hold the +U potential at the occupations it has NOW, through every later refresh (\a frozen), or
    //! let it follow the density again.  Linear-response U is defined frozen (Timrov eq 20); so is the
    //! finite-difference cross-check at \f$U_{\rm in}\ne0\f$.  A frozen term's response kernel is zero.
    virtual void FreezeOccupations(bool frozen) = 0;
    virtual bool OccupationsFrozen() const = 0;
};
} // namespace
