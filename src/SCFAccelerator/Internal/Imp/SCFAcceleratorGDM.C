// File: SCFAcceleratorGDM.C  Implementation of Grassmann-manifold direct minimization.
module;
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <blaze/Math.h>
module qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;

namespace qchem::SCFAccelerators
{
using std::cout;
using std::endl;

static rsmat_t sympart(const rmat_t& B)
{
    size_t n=B.rows();
    rsmat_t S(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<=i;j++)
            S(i,j)=0.5*(B(i,j)+B(j,i));
    return S;
}

// Parallel-transport operator along the Grassmann geodesic whose tangent has SVD
// H = U diag(theta) Vt and start point Y (n x no orthonormal).  Edelman-Arias-Smith eq 2.65:
//   T = (-Y V sin(theta) + U cos(theta)) U^T + (I - U U^T),   V = trans(Vt).
// A tangent W at Y is transported to T*W at the geodesic endpoint.
static rmat_t TransportOp(const rmat_t& Y, const rmat_t& U, const rvec_t& theta, const rmat_t& Vt)
{
    size_t n=Y.rows(), no=Y.columns();
    blaze::DiagonalMatrix<rmat_t> Dc(no),Ds(no);
    for (size_t k=0;k<no;k++){ Dc(k,k)=std::cos(theta[k]); Ds(k,k)=std::sin(theta[k]); }
    rmat_t A = -Y*trans(Vt)*Ds + U*Dc;       // n x no
    rmat_t T = A*trans(U) - U*trans(U);       // n x n
    for (size_t k=0;k<n;k++) T(k,k)+=1.0;     // + I  (the (I - U U^T) part)
    return T;
}

SCFIrrepAcceleratorGDM::SCFIrrepAcceleratorGDM(const GDMParams& p,const LASolver<double>* las,const Irrep& ir,int occ)
: itsParams(p), itsLASolver(las), itsIrrep(ir), itsNocc(occ), itsHaveC(false), itsEn(0.0), itsActive(false)
{
    assert(itsLASolver);
}
SCFIrrepAcceleratorGDM::~SCFIrrepAcceleratorGDM() {}

void SCFIrrepAcceleratorGDM::UseFD(const rsmat_t& F, const rsmat_t& DPrime)
{
    itsFp = itsLASolver->Transform(F);          // F' = Vd F V (orthonormal basis)
    rmat_t E = itsFp*DPrime - DPrime*itsFp;      // [F',D'] commutator (anti-symmetric)
    itsEn = norm(E);
}

LASolver<double>::UUd_t SCFIrrepAcceleratorGDM::NextOrbitals()
{
    if (!ComputeStep())   // seed / diagonalize case
    {
        auto t   = itsLASolver->SolveOrtho(itsFp);
        itsCp    = std::get<1>(t);               // orthonormal-basis orbitals (n x n)
        itsHaveC = true;
        itsHavePrev = false;                     // restart CG after any diagonalizing step
        return t;
    }
    return OrbitalsAt(1.0,true);                 // take the full (diagonal-model) step and commit
}

// Compute the search direction and geodesic for the current Fock, WITHOUT moving the
// orbitals.  Caches everything OrbitalsAt() needs (geodesic SVD, pc orbitals/energies and
// the pending CG history).  Returns false in the seed/diagonalize case (no step computed).
bool SCFIrrepAcceleratorGDM::ComputeStep()
{
    itsActive=false;
    size_t n=itsFp.rows(), no=itsNocc;
    if (!itsHaveC || no==0 || no>=n || itsEn>=itsParams.EMax) return false;
    itsActive=true;
    size_t nv=n-no;

    // Fock in the current MO (orthonormal) basis, and its blocks.
    rmat_t FMO  = trans(itsCp)*itsFp*itsCp;
    rmat_t Cocc = blaze::submatrix(itsCp,0,0 ,n,no);
    rmat_t Cvir = blaze::submatrix(itsCp,0,no,n,nv);
    rsmat_t Foo = sympart(blaze::submatrix(FMO,0 ,0 ,no,no));
    rsmat_t Fvv = sympart(blaze::submatrix(FMO,no,no,nv,nv));
    rmat_t Fov  = blaze::submatrix(FMO,no,0,nv,no);

    // (1) Pseudo-canonicalize: diagonalize occ-occ and virt-virt blocks.
    rvec_t eo, ev; rmat_t Ro, Rv;
    blaze::eigen(Foo, eo, Ro);
    blaze::eigen(Fvv, ev, Rv);
    rmat_t CoccPC = Cocc*Ro;
    rmat_t CvirPC = Cvir*Rv;
    rmat_t g = trans(Rv)*Fov*Ro;        // nv x no  orbital gradient (o-v block)

    // (2) Preconditioned gradient; lift gradient and PG to full ambient tangents at CoccPC.
    rmat_t PG(nv,no);
    for (size_t a=0;a<nv;a++)
        for (size_t i=0;i<no;i++)
            PG(a,i) = g(a,i)/(ev[a]-eo[i]);
    rmat_t Gfull  = CvirPC*g;
    rmat_t PGfull = CvirPC*PG;
    double denom  = sum(g % PG);

    // (3) Conjugate-gradient direction (Polak-Ribiere) with parallel transport.
    rmat_t Dfull;
    if (itsHavePrev)
    {
        rmat_t T   = TransportOp(itsYprev,itsUgeo,itsSgeo,itsVtgeo);
        rmat_t PGt = T*itsPGprev*Ro;
        rmat_t Dt  = T*itsDprev *Ro;
        double beta = (sum(Gfull % PGfull) - sum(Gfull % PGt))/itsDenomPrev;
        if (beta<0.0) beta=0.0;
        Dfull = -PGfull + beta*Dt;
        if (sum(Gfull % Dfull) >= 0.0) Dfull = -PGfull;
    }
    else
        Dfull = -PGfull;

    // (4) Default step length from the diagonal quadratic model:  t* = -<g,d>/<d,H d>.
    rmat_t d = trans(CvirPC)*Dfull;
    double gd=sum(g % d), dHd=0.0;
    for (size_t a=0;a<nv;a++)
        for (size_t i=0;i<no;i++) dHd += d(a,i)*d(a,i)*(ev[a]-eo[i]);
    double t = (dHd>0.0) ? -gd/dHd : 1.0;
    if (t<=0.0) { Dfull=-PGfull; d=trans(CvirPC)*Dfull; t=1.0; }

    // (5) SVD of the tangent H=Dfull -> geodesic factors.
    rmat_t H = CvirPC*d;
    rmat_t U,Vt; rvec_t s;
    blaze::svd(H,U,s,Vt);

    itsScoccPC=CoccPC; itsScvirPC=CvirPC; itsSU=U; itsSVt=Vt; itsSs=s;
    itsSeo=eo; itsSev=ev; itsSPGfull=PGfull; itsSDfull=Dfull; itsSdenom=denom; itsStdef=t;
    return true;
}

// Orbitals at geodesic fraction t of the default step (subject to the trust radius).  With
// commit=true the orbitals and the CG history are updated (a real step); with commit=false
// it is a pure trial used to evaluate the energy in a line search.
LASolver<double>::UUd_t SCFIrrepAcceleratorGDM::OrbitalsAt(double t, bool commit)
{
    size_t n=itsCp.rows(), no=itsNocc, nv=n-no;
    rvec_t s=itsSs;
    double frac = itsStdef*t;          // t is a fraction of the default (quadratic-model) step
    // Trust radius: cap the largest principal rotation angle (cf. ComputeStep).
    double amax = (s.size()>0) ? max(s)*frac : 0.0;
    if (amax>itsParams.Trust) frac *= itsParams.Trust/amax;
    rvec_t theta = s*frac;
    blaze::DiagonalMatrix<rmat_t> Dc(no),Ds(no);
    for (size_t k=0;k<no;k++){ Dc(k,k)=std::cos(theta[k]); Ds(k,k)=std::sin(theta[k]); }
    rmat_t Co = itsScoccPC*trans(itsSVt)*Dc*itsSVt + itsSU*Ds*itsSVt;   // new occupied

    // Complete to a full orthonormal set: virtuals orthogonal to the new occupied block.
    rmat_t W0 = itsScvirPC - Co*(trans(Co)*itsScvirPC);
    rmat_t Uw,Vtw; rvec_t sw;
    blaze::svd(W0,Uw,sw,Vtw);

    rmat_t Cnew(n,n);
    blaze::submatrix(Cnew,0,0 ,n,no) = Co;
    blaze::submatrix(Cnew,0,no,n,nv) = Uw;
    rvec_t e(n);
    for (size_t i=0;i<no;i++) e[i]    = itsSeo[i];
    for (size_t a=0;a<nv;a++) e[no+a] = itsSev[a];

    if (commit)
    {
        itsPGprev=itsSPGfull; itsDprev=itsSDfull; itsDenomPrev=itsSdenom; itsYprev=itsScoccPC;
        itsUgeo=itsSU; itsVtgeo=itsSVt; itsSgeo=theta; itsHavePrev=true;
        itsCp = Cnew;
    }
    rmat_t Uao = itsLASolver->BackTransform(Cnew);
    return std::make_tuple(Uao,Cnew,e);
}

//----------------------------------------------------------------------------------------
// Non-irrep (aggregate) code
//
SCFAcceleratorGDM::SCFAcceleratorGDM(const GDMParams& p) : itsParams(p), itsEn(0.0) {}
SCFAcceleratorGDM::~SCFAcceleratorGDM() {}

SCFIrrepAccelerator* SCFAcceleratorGDM::Create(const LASolver<double>* las,const Irrep& qns,int occ)
{
    if (occ>0)
    {
        itsIrreps.push_back(new SCFIrrepAcceleratorGDM(itsParams,las,qns,occ));
        return itsIrreps.back();
    }
    return new SCFIrrepAcceleratorNull(las,qns);
}

double SCFAcceleratorGDM::GetError() const
{
    double e=0.0;
    for (auto k:itsIrreps) e=std::max(e,k->GetError());
    return e;
}

void SCFAcceleratorGDM::ShowLabels(std::ostream& os) const { os << "  |∇|  Nactive"; }
void SCFAcceleratorGDM::ShowConvergence(std::ostream& os) const
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

} //namespace
