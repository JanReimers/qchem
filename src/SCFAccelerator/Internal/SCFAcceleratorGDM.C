// File: SCFAcceleratorGDM.C  Geometric Direct Minimization (Grassmann-manifold SCF).
//
// First cut: preconditioned steepest descent on the Grassmann manifold, per irrep.
//   - gradient            g   = occupied-virtual block of the (pseudo-canonical) Fock matrix
//   - preconditioned step  d   = -g_ai / (eps_a - eps_i)   (diagonal orbital Hessian)
//   - geodesic step (t=1)  via the SVD of the tangent (Edelman-Arias-Smith geodesic)
// The preconditioner is the inverse diagonal Hessian, so the natural step length is t=1
// (a Newton step), which sidesteps a line search.  Conjugate-gradient + parallel
// transport can be layered on top of this later for faster convergence.
//
// TEMPLATED on T (double | dcmplx): the whole Grassmann machinery generalises to the COMPLEX
// Stiefel/Grassmann manifold by replacing every transpose with a conjugate-transpose and every
// tangent-space inner product <A,B> with its real part Re Tr(A^H B).  For T=double conj is a
// no-op, so tSCFAcceleratorGDM<double> is bit-identical to the original real code.  The complex
// instantiation is what lets a plane-wave / GPW SCF at a genuinely complex k polish past a DIIS
// limit cycle (and converge a degenerate open shell), which DIIS alone cannot.
//
// Ref: Van Voorhis & Head-Gordon, Mol. Phys. 100, 1713 (2002);
//      Edelman, Arias & Smith, SIAM J. Matrix Anal. Appl. 20, 303 (1998).
module;
#include <iosfwd>
#include <vector>
export module qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
export import qchem.SCFAccelerator;
import qchem.Blaze;

export namespace qchem::SCFAccelerators
{

// (GDMParams moved to the public qchem.SCFAccelerator 2026-08-09 -- documented there.)


//! \brief The SCALAR AGGREGATION face one per-irrep GDM block shows its manager (doc/RealComplexPlan.md
//! §6): the manager only ever folds scalar answers (error norms, readiness, engageability, step
//! rejection) over its blocks, so it needs no view of the block scalar T -- one non-template manager
//! serves real and complex blocks alike.
class GDM_Block
{
public:
    virtual ~GDM_Block() {};
    virtual double GetError()   const=0;   //!< this block's ||[F',D']||
    virtual bool   Ready()      const=0;   //!< geodesic step ready THIS iteration (seeded, gated, not forced-diag)
    virtual bool   Engageable() const=0;   //!< STANDING precondition (idempotent + leading-block D')
    virtual bool   Active()     const=0;   //!< contributing to the current step (the Nactive display)
    virtual bool   RejectStep()      =0;   //!< bail-out (same final overrider as tSCFIrrepAccelerator's)
};

template <class T> class tSCFIrrepAcceleratorGDM
    : public virtual tSCFIrrepAccelerator<T>
    , public virtual GDM_Block
{
public:
    tSCFIrrepAcceleratorGDM(const GDMParams&,const LASolver<T>*,const Irrep&,int occ);
    virtual ~tSCFIrrepAcceleratorGDM();

    virtual void UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime);
    virtual typename LASolver<T>::UUd_t NextOrbitals();
    // Direct-minimization line-search hooks (see the base interface):
    //   ComputeStep() does gradient -> CG direction -> geodesic SVD (no orbital change), and
    //     returns false in the seed/diagonalize case (the caller should diagonalize);
    //   OrbitalsAt(t,commit) returns the orbitals at geodesic fraction t of the default step.
    // NextOrbitals() == ComputeStep() ? OrbitalsAt(default_t,true) : diagonalize.
    virtual bool ComputeStep();
    virtual typename LASolver<T>::UUd_t OrbitalsAt(double t, bool commit);
    //! Bail-out callback (see the base contract): drop the CG history, shrink the trust radius, rebuild a
    //! steepest-descent direction, and report whether a retry is worth it.  false => EXHAUSTED, and this
    //! object has armed itsForceDiag so ComputeStep()/Ready() go false and the caller's fallback is safe.
    virtual bool RejectStep();
private:
    // The GDM_Block scalar face (private overrides; the manager calls through GDM_Block*).
    virtual double GetError() const {return itsEn;}
    //! Ready to take a geodesic step THIS iteration (mirrors ComputeStep's gate, minus the Fock build):
    //! seeded, a well-posed occ/virt split, and the residual [F,D] already below FDMax.  itsFp is empty
    //! before the first UseFD, so this is false (unsigned itsNocc < 0 rows) until seeded -- as intended.
    //! (itsForceDiag: armed by an EXHAUSTED RejectStep, so both this and ComputeStep() report "no geodesic"
    //! until a diagonalizing step has actually been taken -- the safety guarantee the bail-out contract owes
    //! the caller.  Cleared in NextOrbitals' diagonalize branch, NOT in ComputeStep, so it cannot be consumed
    //! by a mere query.)
    virtual bool Ready() const { return itsHaveC && itsNocc>0 && itsNocc<itsFp.rows() && itsEn<itsParams.FDMax
                                && !itsForceDiag && itsBlockOccupied && itsIdempotent; }
    //! The STANDING precondition (see tSCFAccelerator::Engageable): an integer-occupation density
    //! (D' idempotent -- false under Fermi smearing, knowable from the very first UseFD with no orbitals)
    //! that occupies THIS minimizer's leading block (false under MOM; knowable only once seeded, so it
    //! stays true -- optimistic -- until then and the ladder's retreat covers the late discovery).
    virtual bool Engageable() const { return itsIdempotent && itsBlockOccupied; }
    virtual bool Active()     const { return itsActive; }

    GDMParams               itsParams;
    const LASolver<T>*      itsLASolver;
    Irrep                   itsIrrep;
    size_t                  itsNocc;
    bool                    itsHaveC;   //Have we cached a set of orbitals yet?
    mat_t<T>                itsCp;      //Orthonormal-basis orbitals (n x n), columns = MOs.
    hmat_t<T>               itsFp;      //Orthonormal-basis Fock matrix.
    double                  itsEn;      //||[F',D']|| error for this irrep.
    bool                    itsActive;
    double                  itsTrust;      //LIVE trust radius (starts at itsParams.Trust; RejectStep shrinks
                                           //  it, an ACCEPTED step restores it -- the backoff is per-step
                                           //  recovery, not a permanent ratchet).
    bool                    itsForceDiag=false; //EXHAUSTED: force a diagonalizing step (see Ready()).
    //! GDM's UNSTATED PRECONDITION, now checked (see UseFD).  Every block below -- Cocc, Cvir, the o-v
    //! gradient, the geodesic -- is carved out of itsCp BY INDEX, i.e. GDM assumes the occupied states are
    //! the LEADING itsNocc columns.  Aufbau makes that true.  MOM (which occupies by overlap, deliberately
    //! NOT the lowest) and Fermi smearing (fractional over all states) make it FALSE, and the gradient is
    //! then built for a different manifold than the density actually occupies -- silently, since nothing
    //! here ever looked at D' beyond its commutator norm.  MEASURED on MnO: the leading-block occupation is
    //! 14.5 Ha ABOVE the MOM/smeared one at the SAME orbitals, so the minimizer was optimising a state the
    //! run had already rejected.  When this is false, GDM declines to engage at all.
    bool                    itsBlockOccupied=true;
    //! Is D' IDEMPOTENT (D'^2 == D', i.e. integer occupations)?  The half of the standing precondition that
    //! needs NO orbitals: computed on every UseFD from D' alone, so the ladder can veto a hand-off BEFORE
    //! this rung ever diagonalizes.  Fermi smearing makes it false (fractional occupations over all states);
    //! a MOM occupation is a genuine projector and passes (itsBlockOccupied is what catches MOM).
    bool                    itsIdempotent=true;

    // Conjugate-gradient state (Polak-Ribiere + parallel transport on the manifold).
    bool     itsHavePrev=false;
    mat_t<T> itsPGprev, itsDprev; //prev preconditioned gradient & search direction (full tangents at itsYprev)
    mat_t<T> itsYprev;            //prev pseudo-canonical occupied orbitals (start of last geodesic)
    mat_t<T> itsUgeo, itsVtgeo;   //SVD factors of last geodesic tangent (H = Ugeo diag Vtgeo)
    rvec_t   itsSgeo;             //last geodesic angles (singular values * step length)
    double   itsDenomPrev=0.0;    //Re<G_old, P G_old> for the Polak-Ribiere denominator

    // Step state cached by ComputeStep() and consumed by OrbitalsAt().
    mat_t<T> itsScoccPC, itsScvirPC; //pseudo-canonical occupied/virtual orbitals
    mat_t<T> itsSU, itsSVt;          //SVD factors of the tangent H = itsSU diag(itsSs) itsSVt
    rvec_t   itsSs, itsSeo, itsSev;  //singular values; pc occupied/virtual orbital energies
    mat_t<T> itsSPGfull, itsSDfull;  //full tangents for the CG history
    double   itsSdenom=0.0;          //Re<G,PG> for this step
    double   itsStdef=1.0;           //default (diagonal quadratic-model) step length
};

//! The GDM MANAGER -- non-template (doc/RealComplexPlan.md §6): "no global coupling, each irrep steps
//! on its own", so its state is a parameter block and a scalar; the typed \c Create overloads hand out
//! block-typed minimizers and register their scalar GDM_Block face.
class SCFAcceleratorGDM : public virtual SCFAccelerator
{
public:
    SCFAcceleratorGDM(const GDMParams&);
    ~SCFAcceleratorGDM();
    virtual tSCFIrrepAccelerator<double>* Create(const LASolver<double>*,const Irrep&, int occ);
    virtual tSCFIrrepAccelerator<dcmplx>* Create(const LASolver<dcmplx>*,const Irrep&, int occ);
    virtual bool   CalculateProjections() {return true;} //No global coupling: each irrep steps on its own.
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;
    virtual const char* Tag() const {return "GDM";}
    virtual bool   WantsLineSearch() const {return true;} //GDM is run by the direct-min loop.
    virtual bool   CanLineSearch() const; //true once every irrep is Ready() (seeded + [F,D]<FDMax)
    virtual bool   Engageable() const;    //every irrep's STANDING precondition (idempotent + leading-block D')
    //! Bail-out: reject the step in EVERY irrep.  A retry is offered only if EVERY irrep can still retry --
    //! the line search takes ONE t across all of them, so a step is only meaningful if they all still have a
    //! usable direction.  One exhausted irrep exhausts the whole accelerator (and each irrep has armed its
    //! own forced diagonalize, so the caller's fallback is safe for all of them).
    virtual bool   RejectStep();
private:
    template <class T> tSCFIrrepAccelerator<T>* CreateT(const LASolver<T>*,const Irrep&, int occ);

    GDMParams itsParams;
    std::vector<GDM_Block*> itsIrreps;   // the scalar aggregation face -- children may differ in T
    double itsEn;
};

using SCFIrrepAcceleratorGDM  = tSCFIrrepAcceleratorGDM<double>;  using cSCFIrrepAcceleratorGDM = tSCFIrrepAcceleratorGDM<dcmplx>;

} //namespace
