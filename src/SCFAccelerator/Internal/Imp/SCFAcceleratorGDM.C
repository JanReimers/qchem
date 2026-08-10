// File: SCFAcceleratorGDM.C  Implementation of Grassmann-manifold direct minimization.
//
// TEMPLATED on T (double | dcmplx).  The COMPLEX generalisation of every step is mechanical but
// must be got exactly right (a stray transpose-instead-of-conjugate-transpose is a no-op for real
// and silently wrong for complex): (a) every M^T becomes M^H (blazem::ctrans); (b) every tangent-
// space inner product <A,B> is the REAL part of the Hermitian Frobenius product Re Tr(A^H B)
// (ReDot below); (c) the diagonal-Hessian denominator sums |d|^2 (Abs2), not d^2.  For T=double
// conj is the identity, so this is bit-identical to the original real GDM.
module;
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
module qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
import qchem.Blaze;

namespace qchem::SCFAccelerators
{
using std::cout;
using std::endl;

// Scalar conjugate that is the IDENTITY for real T (std::conj(double) would promote to complex).
static inline double Conj(double x){return x;}
static inline dcmplx Conj(dcmplx z){return std::conj(z);}
// Squared magnitude |x|^2 as a real, for both T.
static inline double Abs2(double x){return x*x;}
static inline double Abs2(dcmplx z){return std::norm(z);}

// Tangent-space inner product Re Tr(A^H B) = Re sum_ij conj(A_ij) B_ij (real for both T; the
// Grassmann metric).  For T=double conj is a no-op so this is the ordinary sum(A % B).
template <class T> static double ReDot(const mat_t<T>& A, const mat_t<T>& B)
{
    return std::real(blazem::sum(blazem::conj(A) % B));
}

// Hermitian part 1/2 (B + B^H) as an hmat_t<T> (symmetric part for real T).  Filled on the lower
// triangle so the HermitianMatrix auto-mirrors the upper as its conjugate and the diagonal is real.
template <class T> static hmat_t<T> hermpart(const mat_t<T>& B)
{
    size_t n=B.rows();
    hmat_t<T> S(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<=i;j++)
            S(i,j)=0.5*(B(i,j)+Conj(B(j,i)));   // diagonal (j==i) -> real by construction
    return S;
}

// Parallel-transport operator along the Grassmann geodesic whose tangent has SVD H = U diag(theta) V^H
// and start point Y (n x no orthonormal).  Edelman-Arias-Smith eq 2.65 (complex form, V = Vt^H):
//   T = (-Y V sin(theta) + U cos(theta)) U^H + (I - U U^H).
// A tangent W at Y is transported to T*W at the geodesic endpoint.
template <class T> static mat_t<T> TransportOp(const mat_t<T>& Y, const mat_t<T>& U, const rvec_t& theta, const mat_t<T>& Vt)
{
    size_t n=Y.rows(), no=Y.columns();
    blazem::DiagonalMatrix<mat_t<T>> Dc(no),Ds(no);
    for (size_t k=0;k<no;k++){ Dc(k,k)=std::cos(theta[k]); Ds(k,k)=std::sin(theta[k]); }
    mat_t<T> A = -Y*blazem::ctrans(Vt)*Ds + U*Dc;               // n x no  (V = Vt^H)
    mat_t<T> Tm = A*blazem::ctrans(U) - U*blazem::ctrans(U);    // n x n
    for (size_t k=0;k<n;k++) Tm(k,k)+=T(1.0);                   // + I  (the (I - U U^H) part)
    return Tm;
}

template <class T> tSCFIrrepAcceleratorGDM<T>::tSCFIrrepAcceleratorGDM(const GDMParams& p,const LASolver<T>* las,const Irrep& ir,int occ)
: itsParams(p), itsLASolver(las), itsIrrep(ir), itsNocc(occ), itsHaveC(false), itsEn(0.0), itsActive(false)
, itsTrust(p.Trust)
{
    assert(itsLASolver);
}
template <class T> tSCFIrrepAcceleratorGDM<T>::~tSCFIrrepAcceleratorGDM() {}

template <class T> void tSCFIrrepAcceleratorGDM<T>::UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime)
{
    itsFp = itsLASolver->Transform(F);           // F' = Vd F V (orthonormal basis)
    mat_t<T> E = itsFp*DPrime - DPrime*itsFp;    // [F',D'] commutator (anti-Hermitian)
    itsEn = std::real(blazem::norm(E));
    // PRECONDITION CHECK (itsBlockOccupied): is the density actually the projector onto our LEADING itsNocc
    // columns?  Tr(D' P_block) counts how many occupied electrons live in the block GDM is about to rotate:
    // exactly itsNocc under aufbau, short of it under MOM (a different set) or Fermi smearing (fractional
    // over all states).  Checked HERE because this is the only place the true D' is in hand.  Cheap beside
    // the commutator above: nocc x O(n^2).
    if (itsHaveC && itsNocc>0 && itsNocc<=itsCp.columns())
    {
        // Tr(D' P_block) = sum_k <c_k| D' |c_k> over the leading itsNocc columns.  Written as an explicit
        // double loop rather than a Blaze expression: the operand is a column VIEW of itsCp and a symmetric/
        // Hermitian matrix product, and spelling out the contraction keeps the conj() convention unambiguous
        // for both the real and the complex instantiation.
        double inBlock=0.0;
        for (size_t k=0;k<itsNocc;k++)
        {
            T acc{};
            for (size_t i=0;i<itsCp.rows();i++)
                for (size_t j=0;j<itsCp.rows();j++)
                    acc += Conj(itsCp(i,k))*DPrime(i,j)*itsCp(j,k);
            inBlock += std::real(acc);
        }
        const bool ok = std::fabs(inBlock-(double)itsNocc) < 1e-6*(double)itsNocc;
        if (itsBlockOccupied && !ok)
            std::cerr << "[GDM] DECLINING to engage: the occupied density is NOT this minimizer's leading "
                      << itsNocc << " orbitals (Tr(D'P_block)=" << inBlock << ").  GDM rotates a leading-index"
                         " block, so a MOM or Fermi-smeared occupation puts its gradient on a different"
                         " manifold than the density.  Falling back to diagonalizing steps." << std::endl;
        itsBlockOccupied = ok;
    }
}

template <class T> typename LASolver<T>::UUd_t tSCFIrrepAcceleratorGDM<T>::NextOrbitals()
{
    if (!ComputeStep())   // seed / diagonalize case
    {
        auto t   = itsLASolver->SolveOrtho(itsFp);
        itsCp    = std::get<1>(t);               // orthonormal-basis orbitals (n x n)
        itsHaveC = true;
        itsHavePrev = false;                     // restart CG after any diagonalizing step
        itsForceDiag = false;                    // the forced diagonalize has now HAPPENED: re-arm the geodesic
        itsTrust = itsParams.Trust;              // and give the fresh direction the full trust radius back
        return t;
    }
    return OrbitalsAt(1.0,true);                 // take the full (diagonal-model) step and commit
}

// Compute the search direction and geodesic for the current Fock, WITHOUT moving the
// orbitals.  Caches everything OrbitalsAt() needs (geodesic SVD, pc orbitals/energies and
// the pending CG history).  Returns false in the seed/diagonalize case (no step computed).
template <class T> bool tSCFIrrepAcceleratorGDM<T>::ComputeStep()
{
    itsActive=false;
    size_t n=itsFp.rows(), no=itsNocc;
    // itsForceDiag: an EXHAUSTED RejectStep armed a diagonalizing step.  Report "no geodesic" WITHOUT
    // clearing the flag -- it is cleared only when a diagonalize has actually happened (NextOrbitals),
    // so a query cannot consume the guarantee the bail-out contract owes the caller.
    if (itsForceDiag) return false;
    // The leading-block precondition (see UseFD): refuse rather than optimise the wrong manifold.  This must
    // gate ComputeStep and not merely Ready(), because NextOrbitals() TAKES the step whenever ComputeStep
    // succeeds -- a check only in Ready() would still let the fallback path commit a full geodesic step.
    if (!itsBlockOccupied) return false;
    if (!itsHaveC || no==0 || no>=n || itsEn>=itsParams.FDMax) return false;
    itsActive=true;
    size_t nv=n-no;

    // Fock in the current MO (orthonormal) basis, and its blocks.
    mat_t<T> FMO  = blazem::ctrans(itsCp)*itsFp*itsCp;
    mat_t<T> Cocc = blazem::submatrix(itsCp,0,0 ,n,no);
    mat_t<T> Cvir = blazem::submatrix(itsCp,0,no,n,nv);
    hmat_t<T> Foo = hermpart<T>(blazem::submatrix(FMO,0 ,0 ,no,no));
    hmat_t<T> Fvv = hermpart<T>(blazem::submatrix(FMO,no,no,nv,nv));
    mat_t<T> Fov  = blazem::submatrix(FMO,no,0,nv,no);

    // (1) Pseudo-canonicalize: diagonalize occ-occ and virt-virt blocks (Hermitian -> real energies).
    rvec_t eo, ev; mat_t<T> Ro, Rv;
    blazem::eigen(Foo, eo, Ro);
    blazem::eigen(Fvv, ev, Rv);
    mat_t<T> CoccPC = Cocc*Ro;
    mat_t<T> CvirPC = Cvir*Rv;
    mat_t<T> g = blazem::ctrans(Rv)*Fov*Ro;      // nv x no  orbital gradient (o-v block)

    // (2) Preconditioned gradient; lift gradient and PG to full ambient tangents at CoccPC.
    // The gap is FLOORED to keep the preconditioner positive-definite (GDMParams::PCFloor): it vanishes on a
    // near-degenerate pair and goes NEGATIVE on a non-aufbau one (a hole), and an unfloored quotient then
    // either explodes or flips sign into an ascent direction.  A preconditioner is a metric, not the Hessian,
    // so flooring rescales stiff modes and never costs the descent property.
    mat_t<T> PG(nv,no);
    for (size_t a=0;a<nv;a++)
        for (size_t i=0;i<no;i++)
            PG(a,i) = g(a,i)/std::max(ev[a]-eo[i], itsParams.PCFloor);
    mat_t<T> Gfull  = CvirPC*g;
    mat_t<T> PGfull = CvirPC*PG;
    double denom  = ReDot<T>(g,PG);

    // (3) Conjugate-gradient direction (Polak-Ribiere) with parallel transport.
    mat_t<T> Dfull;
    if (itsHavePrev)
    {
        mat_t<T> Tr  = TransportOp<T>(itsYprev,itsUgeo,itsSgeo,itsVtgeo);
        mat_t<T> PGt = Tr*itsPGprev*Ro;
        mat_t<T> Dt  = Tr*itsDprev *Ro;
        double beta = (ReDot<T>(Gfull,PGfull) - ReDot<T>(Gfull,PGt))/itsDenomPrev;
        if (beta<0.0) beta=0.0;
        Dfull = -PGfull + beta*Dt;
        if (ReDot<T>(Gfull,Dfull) >= 0.0) Dfull = -PGfull;
    }
    else
        Dfull = -PGfull;

    // (4) Default step length from the diagonal quadratic model:  t* = -<g,d>/<d,H d>.
    mat_t<T> d = blazem::ctrans(CvirPC)*Dfull;
    double gd=ReDot<T>(g,d), dHd=0.0;
    // Same floor as the preconditioner above: an indefinite <d,Hd> makes the quadratic-model step length
    // -gd/dHd meaningless (the dHd>0 guard below then discards it entirely and falls back to t=1), so a
    // hole would cost the model step rather than merely bias it.
    for (size_t a=0;a<nv;a++)
        for (size_t i=0;i<no;i++) dHd += Abs2(d(a,i))*std::max(ev[a]-eo[i], itsParams.PCFloor);
    double t = (dHd>0.0) ? -gd/dHd : 1.0;
    if (t<=0.0) { Dfull=-PGfull; d=blazem::ctrans(CvirPC)*Dfull; t=1.0; }

    // (5) SVD of the tangent H=Dfull -> geodesic factors.
    mat_t<T> H = CvirPC*d;
    mat_t<T> U,Vt; rvec_t s;
    blazem::svd(H,U,s,Vt);

    itsScoccPC=CoccPC; itsScvirPC=CvirPC; itsSU=U; itsSVt=Vt; itsSs=s;
    itsSeo=eo; itsSev=ev; itsSPGfull=PGfull; itsSDfull=Dfull; itsSdenom=denom; itsStdef=t;
    return true;
}

// BAIL-OUT (the RejectStep contract): the caller evaluated the step we proposed and it was NOT acceptable.
// Two things are wrong with a rejected direction and we address both.  (1) The CONJUGACY HISTORY: a CG
// direction is built from the previous one, so one bad direction poisons every direction after it -- drop it
// and fall back to steepest descent, which is the direction a line search can always make progress along IF
// the energy is a continuous function of t.  (2) The STEP SCALE: shrink the trust radius, so the retry probes
// a genuinely smaller neighbourhood instead of re-proposing the same rotation the caller just rejected.
//
// Backing off has a floor.  Below TrustMin we are no longer testing the direction, we are testing whether
// E(t) -> E(0) at all -- and it does not when the CALLER's occupation decision jumps branch between trial
// points (measured on MnO: the best of 12 backtracks still +14.5 Ha at t=2.4e-4, which is a DISCONTINUITY in
// t, not a bad direction; doc/SymmetryUpgradePlan.md 7 step 7).  No trust radius can fix that, so we stop
// pretending and declare EXHAUSTED -- arming itsForceDiag, which is what makes the caller's fallback safe.
template <class T> bool tSCFIrrepAcceleratorGDM<T>::RejectStep()
{
    itsHavePrev = false;                 // (1) the CG history produced a rejected direction -- restart it
    itsTrust   *= itsParams.TrustBackoff;// (2) and probe a smaller neighbourhood
    if (itsTrust < itsParams.TrustMin)
    {
        itsForceDiag = true;             // EXHAUSTED: Ready()/ComputeStep() now false => caller may fall back
        return false;
    }
    return ComputeStep();                // rebuild (steepest descent, tighter trust); false => also exhausted
}

// Orbitals at geodesic fraction t of the default step (subject to the trust radius).  With
// commit=true the orbitals and the CG history are updated (a real step); with commit=false
// it is a pure trial used to evaluate the energy in a line search.
template <class T> typename LASolver<T>::UUd_t tSCFIrrepAcceleratorGDM<T>::OrbitalsAt(double t, bool commit)
{
    size_t n=itsCp.rows(), no=itsNocc, nv=n-no;
    rvec_t s=itsSs;
    double frac = itsStdef*t;          // t is a fraction of the default (quadratic-model) step
    // Trust radius: cap the largest principal rotation angle (cf. ComputeStep).
    double amax = (s.size()>0) ? blazem::max(s)*frac : 0.0;
    if (amax>itsTrust) frac *= itsTrust/amax;
    rvec_t theta = s*frac;
    blazem::DiagonalMatrix<mat_t<T>> Dc(no),Ds(no);
    for (size_t k=0;k<no;k++){ Dc(k,k)=std::cos(theta[k]); Ds(k,k)=std::sin(theta[k]); }
    // Geodesic Y(t) = Y V cos(theta) V^H + U sin(theta) V^H, with V = itsSVt^H (SVD H = U s itsSVt).
    mat_t<T> Co = itsScoccPC*blazem::ctrans(itsSVt)*Dc*itsSVt + itsSU*Ds*itsSVt;   // new occupied

    // Complete to a full orthonormal set: virtuals orthogonal to the new occupied block.
    mat_t<T> W0 = itsScvirPC - Co*(blazem::ctrans(Co)*itsScvirPC);
    mat_t<T> Uw,Vtw; rvec_t sw;
    blazem::svd(W0,Uw,sw,Vtw);

    mat_t<T> Cnew(n,n);
    blazem::submatrix(Cnew,0,0 ,n,no) = Co;
    blazem::submatrix(Cnew,0,no,n,nv) = Uw;
    rvec_t e(n);
    for (size_t i=0;i<no;i++) e[i]    = itsSeo[i];
    for (size_t a=0;a<nv;a++) e[no+a] = itsSev[a];

    if (commit)
    {
        itsTrust = itsParams.Trust;   // an ACCEPTED step earns the full trust radius back (see RejectStep)
        itsPGprev=itsSPGfull; itsDprev=itsSDfull; itsDenomPrev=itsSdenom; itsYprev=itsScoccPC;
        itsUgeo=itsSU; itsVtgeo=itsSVt; itsSgeo=theta; itsHavePrev=true;
        itsCp = Cnew;
    }
    mat_t<T> Uao = itsLASolver->BackTransform(Cnew);
    return std::make_tuple(Uao,Cnew,e);
}

//----------------------------------------------------------------------------------------
// Non-irrep (aggregate) code
//
template <class T> tSCFAcceleratorGDM<T>::tSCFAcceleratorGDM(const GDMParams& p) : itsParams(p), itsEn(0.0) {}
template <class T> tSCFAcceleratorGDM<T>::~tSCFAcceleratorGDM() {}

template <class T> tSCFIrrepAccelerator<T>* tSCFAcceleratorGDM<T>::Create(const LASolver<T>* las,const Irrep& qns,int occ)
{
    if (occ>0)
    {
        itsIrreps.push_back(new tSCFIrrepAcceleratorGDM<T>(itsParams,las,qns,occ));
        return itsIrreps.back();
    }
    return new tSCFIrrepAcceleratorNull<T>(las,qns);
}

template <class T> double tSCFAcceleratorGDM<T>::GetError() const
{
    double e=0.0;
    for (auto k:itsIrreps) e=std::max(e,k->GetError());
    return e;
}

// The direct-min driver runs only when EVERY irrep can take a geodesic step; if any is still seeding or
// above FDMax, the iterator falls back to a stable MIXED fixed-point step instead.
template <class T> bool tSCFAcceleratorGDM<T>::RejectStep()
{
    // Reject in EVERY irrep -- they all contributed to the step the caller rejected, and every one of them
    // must arm its own forced diagonalize if it is exhausted, or the fallback is unsafe for that block.
    // Deliberately NOT short-circuited: `all &= k->RejectStep()` would skip the rest on the first false.
    bool all = !itsIrreps.empty();
    for (auto k:itsIrreps) { const bool retry=k->RejectStep(); all = all && retry; }
    // A retry is offered only if EVERY irrep can still retry: the line search takes ONE t across all of
    // them, so a retry is only meaningful when they all still hold a usable direction.
    return all;
}

template <class T> bool tSCFAcceleratorGDM<T>::CanLineSearch() const
{
    if (itsIrreps.empty()) return false;
    for (auto k:itsIrreps) if (!k->Ready()) return false;
    return true;
}

template <class T> void tSCFAcceleratorGDM<T>::ShowLabels(std::ostream& os) const { os << "  |∇|  Nactive"; }
template <class T> void tSCFAcceleratorGDM<T>::ShowConvergence(std::ostream& os) const
{
    os << std::scientific << GetError() << " ";
    size_t nactive=0;
    for (auto k:itsIrreps)
        if(k->itsActive) nactive++;

    if (nactive>0)
        os << std::setw(4) << nactive;
    else
        os << "    ";

    os << "    ";
}

template class tSCFIrrepAcceleratorGDM<double>;
template class tSCFIrrepAcceleratorGDM<dcmplx>;
template class tSCFAcceleratorGDM<double>;
template class tSCFAcceleratorGDM<dcmplx>;

} //namespace
