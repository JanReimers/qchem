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
    itsActive=false;
    size_t n  = itsFp.rows();
    size_t no = itsNocc;
    // First step (or trivial irreps, or still far from convergence): diagonalize & cache.
    if (!itsHaveC || no==0 || no>=n || itsEn>=itsParams.EMax)
    {
        auto t   = itsLASolver->SolveOrtho(itsFp);
        itsCp    = std::get<1>(t);               // orthonormal-basis orbitals (n x n)
        itsHaveC = true;
        return t;
    }
    itsActive=true;
    size_t nv = n - no;

    // Fock in the current MO (orthonormal) basis, and its blocks.
    rmat_t FMO  = trans(itsCp)*itsFp*itsCp;                 // n x n
    rmat_t Cocc = blaze::submatrix(itsCp,0,0 ,n,no);        // n x no
    rmat_t Cvir = blaze::submatrix(itsCp,0,no,n,nv);        // n x nv
    rsmat_t Foo = sympart(blaze::submatrix(FMO,0 ,0 ,no,no));
    rsmat_t Fvv = sympart(blaze::submatrix(FMO,no,no,nv,nv));
    rmat_t Fov = blaze::submatrix(FMO,no,0,nv,no);          // nv x no

    // (1) Pseudo-canonicalize: diagonalize occ-occ and virt-virt blocks.
    rvec_t eo, ev; rmat_t Ro, Rv;
    blaze::eigen(Foo, eo, Ro);                 // Foo = Ro diag(eo) Ro^T (eo ascending)
    blaze::eigen(Fvv, ev, Rv);
    rmat_t CoccPC = Cocc*Ro;            // n x no  pseudo-canonical occupied
    rmat_t CvirPC = Cvir*Rv;            // n x nv  pseudo-canonical virtual
    rmat_t g = trans(Rv)*Fov*Ro;        // nv x no  orbital gradient (o-v block)

    // (2) Preconditioned step d = -g/(ev-eo); the preconditioner is the inverse diagonal
    //     Hessian, so the natural (Newton) step length is t=1 -- no line search needed.
    rmat_t d(nv,no);
    for (size_t a=0;a<nv;a++)
        for (size_t i=0;i<no;i++)
            d(a,i) = -g(a,i)/(ev[a]-eo[i]);

    // (3) Grassmann geodesic at t=1 via the SVD of the tangent H = CvirPC d.
    rmat_t H = CvirPC*d;                 // n x no  tangent (in the virtual space)
    rmat_t U,Vt; rvec_t s;
    blaze::svd(H,U,s,Vt);                       // H = U diag(s) Vt ; U:n x no, Vt:no x no
    blaze::DiagonalMatrix<rmat_t> Dc(no),Ds(no);
    for (size_t k=0;k<no;k++){ Dc(k,k)=std::cos(s[k]); Ds(k,k)=std::sin(s[k]); }
    rmat_t Co = CoccPC*trans(Vt)*Dc*Vt + U*Ds*Vt;          // n x no  new occupied

    // Complete to a full orthonormal set: virtuals orthogonal to the new occupied block.
    rmat_t W0 = CvirPC - Co*(trans(Co)*CvirPC);            // n x nv
    rmat_t Uw,Vtw; rvec_t sw;
    blaze::svd(W0,Uw,sw,Vtw);                                     // Uw: n x nv orthonormal

    rmat_t Cnew(n,n);
    blaze::submatrix(Cnew,0,0 ,n,no) = Co;
    blaze::submatrix(Cnew,0,no,n,nv) = Uw;
    rvec_t e(n);
    for (size_t i=0;i<no;i++) e[i]    = eo[i];
    for (size_t a=0;a<nv;a++) e[no+a] = ev[a];

    itsCp = Cnew;
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
