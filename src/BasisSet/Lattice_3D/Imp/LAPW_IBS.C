// File: BasisSet/Lattice_3D/Imp/LAPW_IBS.C  LAPW Hamiltonian/overlap/nuclear assembly.
module;
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.LAPW_IBS;
import qchem.Symmetry.Factory;   // BlochFactory
import qchem.Math;               // Pi, FourPi, Cube
import qchem.SpecialFunctions;
import qchem.Blaze;              // zeroH
import qchem.Vector3D;           // dot product (operator*)

namespace
{
// Radial-function table on the grid r[0..NQ] (r[0]=0, r[NQ]=R): the muffin-tin radial function
// u_l(r)=R_l(r) at the linearization energy, its r-derivative u_l', the energy derivative udot_l, and
// udot_l' -- plus their boundary values at R.  V(r)=-Z/r inside the sphere (Z=0 => free, analytic).
struct RadTab { rvec_t u,up,ud,udp; double uR,upR,udR,udpR; };

// Regular radial solution u_l=R_l and u_l' on the grid (V=-Z/r, energy E).  Z=0: analytic Bessel
// j_l(qr).  Z!=0: RK4 integration of P''=f P with P=r R_l, f(r)=l(l+1)/r^2 - 2Z/r - 2E, started from
// the r^{l+1} series at r[1] (regular at the origin).
void SolveRadial(int l,double E,double Z,const rvec_t& r,double R,
                 rvec_t& u,rvec_t& up,double& uR,double& upR)
{
    int NQ=static_cast<int>(r.size())-1;
    u=rvec_t(NQ+1,0.0); up=rvec_t(NQ+1,0.0);
    if (Z==0.0)
    {
        double q=std::sqrt(2.0*E);
        for (int i=1;i<=NQ;i++){ double x=q*r[i];
            rvec_t j=SpecialFunctions::SphericalBessel(l,x);
            rvec_t jp=SpecialFunctions::SphericalBesselPrime(l,x,j);
            u[i]=j[l]; up[i]=q*jp[l]; }
        rvec_t jR=SpecialFunctions::SphericalBessel(l,q*R);
        rvec_t jpR=SpecialFunctions::SphericalBesselPrime(l,q*R,jR);
        uR=jR[l]; upR=q*jpR[l];
        return;
    }
    auto f=[&](double rr){ return l*(l+1)/(rr*rr) - 2.0*Z/rr - 2.0*E; };
    double h=R/NQ;
    double P=std::pow(r[1],l+1), Q=(l+1)*std::pow(r[1],l);   // u ~ r^l series start at r[1]
    u[1]=P/r[1]; up[1]=Q/r[1]-P/(r[1]*r[1]);
    for (int i=1;i<NQ;i++)
    {
        double rr=r[i];
        double k1P=Q,           k1Q=f(rr)*P;
        double k2P=Q+0.5*h*k1Q, k2Q=f(rr+0.5*h)*(P+0.5*h*k1P);
        double k3P=Q+0.5*h*k2Q, k3Q=f(rr+0.5*h)*(P+0.5*h*k2P);
        double k4P=Q+h*k3Q,     k4Q=f(rr+h)*(P+h*k3P);
        P+=h/6.0*(k1P+2*k2P+2*k3P+k4P);
        Q+=h/6.0*(k1Q+2*k2Q+2*k3Q+k4Q);
        double rn=r[i+1];
        u[i+1]=P/rn; up[i+1]=Q/rn-P/(rn*rn);
    }
    uR=u[NQ]; upR=up[NQ];
}

// Build u_l, u_l', udot_l, udot_l' (udot = d/dE).  Z=0: analytic.  Z!=0: u from the ODE at E, udot by
// central finite difference of the ODE solution in E.
RadTab BuildRadial(int l,double E,double Z,const rvec_t& r,double R)
{
    int NQ=static_cast<int>(r.size())-1;
    RadTab t;
    if (Z==0.0)
    {
        double q=std::sqrt(2.0*E);
        t.u=rvec_t(NQ+1,0.0); t.up=rvec_t(NQ+1,0.0); t.ud=rvec_t(NQ+1,0.0); t.udp=rvec_t(NQ+1,0.0);
        for (int i=1;i<=NQ;i++){ double rr=r[i],x=q*rr;
            rvec_t j=SpecialFunctions::SphericalBessel(l,x);
            rvec_t jp=SpecialFunctions::SphericalBesselPrime(l,x,j);
            double jlpp=-(2.0/x)*jp[l]-(1.0-l*(l+1)/(x*x))*j[l];
            t.u[i]=j[l]; t.up[i]=q*jp[l]; t.ud[i]=(rr/q)*jp[l]; t.udp[i]=(1.0/q)*jp[l]+rr*jlpp; }
        double xR=q*R;
        rvec_t jR=SpecialFunctions::SphericalBessel(l,xR);
        rvec_t jpR=SpecialFunctions::SphericalBesselPrime(l,xR,jR);
        double jRpp=-(2.0/xR)*jpR[l]-(1.0-l*(l+1)/(xR*xR))*jR[l];
        t.uR=jR[l]; t.upR=q*jpR[l]; t.udR=(R/q)*jpR[l]; t.udpR=(1.0/q)*jpR[l]+R*jRpp;
        return t;
    }
    double d=1e-4;
    rvec_t u0,up0,up_,upp,um,upm; double uR0,upR0,uRp,upRp,uRm,upRm;
    SolveRadial(l,E,    Z,r,R,u0,up0,uR0,upR0);
    SolveRadial(l,E+d,  Z,r,R,up_,upp,uRp,upRp);
    SolveRadial(l,E-d,  Z,r,R,um,upm,uRm,upRm);
    t.u=u0; t.up=up0; t.uR=uR0; t.upR=upR0;
    t.ud=rvec_t(NQ+1,0.0); t.udp=rvec_t(NQ+1,0.0);
    for (int i=0;i<=NQ;i++){ t.ud[i]=(up_[i]-um[i])/(2*d); t.udp[i]=(upp[i]-upm[i])/(2*d); }
    t.udR=(uRp-uRm)/(2*d); t.udpR=(upRp-upRm)/(2*d);
    return t;
}
} // anon namespace

namespace BasisSet::Lattice_3D
{

LAPW_IBS::LAPW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
                   double Ecut, double Rmt, size_t lmax, double Elin, double Z)
    : BasisSet::IrrepBasisSetImp<dcmplx>(Symmetry::BlochFactory(N,kIndex))
    , itsRecip(recip)
    , itsk(kIndex.x/static_cast<double>(N.x),
           kIndex.y/static_cast<double>(N.y),
           kIndex.z/static_cast<double>(N.z))
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
    , itsRmt(Rmt)
    , itsLmax(lmax)
    , itsElin(Elin)
    , itsZ(Z)
{
    const UnitCell& B=itsRecip.GetCell();
    double Gmax=sqrt(2*Ecut)+B.GetDistance(itsk);
    for (const ivec3_t& m : itsRecip.GetGVectors(Gmax))
    {
        double kG=B.GetDistance(itsk+m);
        if (0.5*kG*kG < Ecut) itsG.push_back(m);
    }
}

void LAPW_IBS::Assemble() const
{
    if (itsCached) return;
    const UnitCell& B=itsRecip.GetCell();
    size_t n=GetNumFunctions();
    int    L=static_cast<int>(itsLmax);
    double R=itsRmt, Omega=itsVolume;
    double Vs=FourPi*R*R*R/3.0;

    // radial grid r[0..NQ]
    const int NQ=800;
    double h=R/NQ;
    rvec_t r(NQ+1);
    for (int i=0;i<=NQ;i++) r[i]=i*h;

    // per-l radial tables, boundary values, and radial integrals (overlap S, kinetic K, nuclear V)
    rvec_t uR(L+1),upR(L+1),udR(L+1),udpR(L+1);
    rvec_t Suu(L+1,0),Suud(L+1,0),Sudud(L+1,0);
    rvec_t Kuu(L+1,0),Kuud(L+1,0),Kudud(L+1,0);
    rvec_t Vuu(L+1,0),Vuud(L+1,0),Vudud(L+1,0);
    for (int l=0;l<=L;l++)
    {
        RadTab t=BuildRadial(l,itsElin,itsZ,r,R);
        uR[l]=t.uR; upR[l]=t.upR; udR[l]=t.udR; udpR[l]=t.udpR;
        for (int i=1;i<=NQ;i++)               // i=0 (r=0) contributes nothing
        {
            double w=(i==NQ ? 1.0 : (i%2 ? 4.0 : 2.0))*h/3.0;   // composite Simpson
            double rr=r[i], r2=rr*rr;
            double u=t.u[i],up=t.up[i],ud=t.ud[i],udp=t.udp[i];
            Suu[l]   += w* u*u*r2;
            Suud[l]  += w* u*ud*r2;
            Sudud[l] += w* ud*ud*r2;
            Kuu[l]   += w*(up*up*r2 + l*(l+1)*u*u);            // <p^2> = int(u'^2 r^2 + l(l+1)u^2)
            Kuud[l]  += w*(up*udp*r2 + l*(l+1)*u*ud);
            Kudud[l] += w*(udp*udp*r2 + l*(l+1)*ud*ud);
            if (itsZ!=0.0)                                     // <V> = int u_a (-Z/r) u_b r^2 dr
            {
                Vuu[l]   += w*(-itsZ)*u*u*rr;
                Vuud[l]  += w*(-itsZ)*u*ud*rr;
                Vudud[l] += w*(-itsZ)*ud*ud*rr;
            }
        }
    }

    // per-plane-wave K=k+G, |K|, and the matching coefficients a_l,b_l (value+derivative at R)
    std::vector<rvec3_t> K(n);
    rvec_t  Kmag(n);
    std::vector<rvec_t> a(n),b(n);
    for (size_t i=0;i<n;i++)
    {
        K[i]=B.ToCartesian(itsk+itsG[i]);
        Kmag[i]=std::sqrt(K[i]*K[i]);
        double KR=Kmag[i]*R;
        rvec_t jK =SpecialFunctions::SphericalBessel(L,KR);
        rvec_t jKp(L+1,0.0);
        if (KR>1e-12) jKp=SpecialFunctions::SphericalBesselPrime(L,KR,jK); // K=0: j_l'(0) contributes 0 below
        a[i]=rvec_t(L+1,0.0); b[i]=rvec_t(L+1,0.0);
        for (int l=0;l<=L;l++)
        {
            double rhs1=jK[l], rhs2=Kmag[i]*jKp[l];   // K=0 -> rhs1=delta_{l0}, rhs2=0
            double W=uR[l]*udpR[l]-upR[l]*udR[l];
            a[i][l]=(rhs1*udpR[l]-rhs2*udR[l])/W;
            b[i][l]=(uR[l]*rhs2-upR[l]*rhs1)/W;
        }
    }

    // assemble <p^2>, O, and <V>  = interstitial + muffin-tin (cached in the mutable members)
    itsKp2 =blazem::zeroH<dcmplx>(n);
    itsO   =blazem::zeroH<dcmplx>(n);
    itsVnuc=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            double OI;
            if (i==j) OI=1.0 - Vs/Omega;
            else { rvec3_t dG=K[i]-K[j]; double dg=std::sqrt(dG*dG);
                   OI=-(FourPi*R*R/Omega)*SpecialFunctions::SphericalBessel1(dg*R)/dg; }
            double KdotK=K[i]*K[j];
            double cosg=(Kmag[i]>1e-12 && Kmag[j]>1e-12) ? KdotK/(Kmag[i]*Kmag[j]) : 1.0;
            cosg=std::max(-1.0,std::min(1.0,cosg));
            rvec_t P=SpecialFunctions::LegendreP(L,cosg);
            double sphereO=0.0, sphereK=0.0, sphereV=0.0;
            for (int l=0;l<=L;l++)
            {
                double ai=a[i][l], bi=b[i][l], aj=a[j][l], bj=b[j][l];
                double pref=(2*l+1)*P[l];
                sphereO += pref*( ai*aj*Suu[l]   + (ai*bj+bi*aj)*Suud[l]   + bi*bj*Sudud[l] );
                sphereK += pref*( ai*aj*Kuu[l]   + (ai*bj+bi*aj)*Kuud[l]   + bi*bj*Kudud[l] );
                sphereV += pref*( ai*aj*Vuu[l]   + (ai*bj+bi*aj)*Vuud[l]   + bi*bj*Vudud[l] );
            }
            itsO   (i,j)=dcmplx(OI + sphereO*FourPi/Omega, 0.0);
            itsKp2 (i,j)=dcmplx(KdotK*OI + sphereK*FourPi/Omega, 0.0); // <p^2> (NO 1/2)
            itsVnuc(i,j)=dcmplx(sphereV*FourPi/Omega, 0.0);            // muffin-tin V; interstitial V=0
        }
    itsCached=true;
}

chmat_t LAPW_IBS::MakeOverlap() const { Assemble(); return itsO; }
chmat_t LAPW_IBS::MakeKinetic() const { Assemble(); return itsKp2; }
chmat_t LAPW_IBS::MakeNuclear(const Structure*) const { Assemble(); return itsVnuc; }

cvec_t LAPW_IBS::operator()(const rvec3_t& r) const
{
    size_t n=GetNumFunctions();
    double inv=1.0/sqrt(itsVolume);
    cvec_t v(n);
    for (size_t i=0; i<n; i++)
    {
        double phase=itsRecip.GetCell().ToCartesian(itsk+itsG[i])*r;
        v[i]=dcmplx(cos(phase),sin(phase))*inv;
    }
    return v;
}

cvec3vec_t LAPW_IBS::Gradient(const rvec3_t& r) const
{
    const dcmplx im(0.0,1.0);
    size_t n=GetNumFunctions();
    double inv=1.0/sqrt(itsVolume);
    cvec3vec_t g(n);
    for (size_t i=0; i<n; i++)
    {
        rvec3_t K=itsRecip.GetCell().ToCartesian(itsk+itsG[i]);
        double phase=K*r;
        dcmplx val=dcmplx(cos(phase),sin(phase))*inv;
        g[i]=vec3_t<dcmplx>(im*K.x*val, im*K.y*val, im*K.z*val);
    }
    return g;
}

std::string LAPW_IBS::BasisSetID() const
{
    return Name()+"|nG="+std::to_string(itsG.size())+"|Rmt="+std::to_string(itsRmt)
                 +"|lmax="+std::to_string(itsLmax)+"|El="+std::to_string(itsElin)
                 +"|Z="+std::to_string(itsZ);
}

std::ostream& LAPW_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " plane waves, Rmt=" << itsRmt
              << " lmax=" << itsLmax << " Elin=" << itsElin << " Z=" << itsZ << ", " << GetSymmetry();
}

} //namespace
