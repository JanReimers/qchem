// File: BasisSet/Lattice_3D/Imp/LAPW_IBS.C  LAPW Hamiltonian/overlap assembly.
module;
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.LAPW_IBS;
import qchem.Symmetry.Factory;   // BlochFactory
import qchem.Math;               // Pi, FourPi, Cube + qchem::Math special functions
import qchem.Blaze;              // zeroH
import qchem.Vector3D;           // dot product (operator*)

namespace BasisSet::Lattice_3D
{

LAPW_IBS::LAPW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
                   double Ecut, double Rmt, size_t lmax, double Elin)
    : BasisSet::IrrepBasisSetImp<dcmplx>(Symmetry::BlochFactory(N,kIndex))
    , itsRecip(recip)
    , itsK(kIndex.x/static_cast<double>(N.x),
           kIndex.y/static_cast<double>(N.y),
           kIndex.z/static_cast<double>(N.z))
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
    , itsRmt(Rmt)
    , itsLmax(lmax)
    , itsElin(Elin)
{
    const UnitCell& B=itsRecip.GetCell();
    double Gmax=sqrt(2*Ecut)+B.GetDistance(itsK);
    for (const ivec3_t& m : itsRecip.GetGVectors(Gmax))
    {
        double kG=B.GetDistance(itsK+m);
        if (0.5*kG*kG < Ecut) itsG.push_back(m);
    }
}

void LAPW_IBS::Assemble(chmat_t& H, chmat_t& O) const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=GetNumFunctions();
    int    L=static_cast<int>(itsLmax);
    double R=itsRmt, Omega=itsVolume;
    double Vs=FourPi*R*R*R/3.0;
    double q=std::sqrt(2.0*itsElin), qR=q*R;

    // --- radial functions and their r-derivatives at the sphere boundary r=R (per l) ---------------
    // u_l(r)=j_l(qr), udot_l(r)=(r/q) j_l'(qr); evaluated/differentiated at R.
    std::vector<double> jq =qchem::Math::SphericalBessel(L,qR);
    std::vector<double> jqp=qchem::Math::SphericalBesselPrime(L,qR,jq);
    std::vector<double> uR(L+1),upR(L+1),udR(L+1),udpR(L+1);
    for (int l=0; l<=L; l++)
    {
        double jlpp=-(2.0/qR)*jqp[l]-(1.0-l*(l+1)/(qR*qR))*jq[l]; // j_l''(qR)
        uR[l]=jq[l];
        upR[l]=q*jqp[l];
        udR[l]=(R/q)*jqp[l];
        udpR[l]=(1.0/q)*jqp[l]+R*jlpp;
    }

    // --- radial integrals over [0,R] by composite Simpson (all integrands vanish at r=0) -----------
    const int NQ=400;                       // even
    double h=R/NQ;
    std::vector<double> Suu(L+1,0.0),Suud(L+1,0.0),Sudud(L+1,0.0);
    std::vector<double> Tuu(L+1,0.0),Tuud(L+1,0.0),Tudud(L+1,0.0);
    for (int iq=1; iq<=NQ; iq++)            // iq=0 (r=0) contributes nothing
    {
        double r=iq*h;
        double w=(iq==NQ ? 1.0 : (iq%2 ? 4.0 : 2.0))*h/3.0;
        double x=q*r;
        std::vector<double> j =qchem::Math::SphericalBessel(L,x);
        std::vector<double> jp=qchem::Math::SphericalBesselPrime(L,x,j);
        double r2=r*r;
        for (int l=0; l<=L; l++)
        {
            double jlpp=-(2.0/x)*jp[l]-(1.0-l*(l+1)/(x*x))*j[l];
            double u =j[l];
            double up=q*jp[l];
            double ud =(r/q)*jp[l];
            double udp=(1.0/q)*jp[l]+r*jlpp;
            Suu[l]   += w* u*u*r2;
            Suud[l]  += w* u*ud*r2;
            Sudud[l] += w* ud*ud*r2;
            Tuu[l]   += w* 0.5*(up*up*r2 + l*(l+1)*u*u);
            Tuud[l]  += w* 0.5*(up*udp*r2 + l*(l+1)*u*ud);
            Tudud[l] += w* 0.5*(udp*udp*r2 + l*(l+1)*ud*ud);
        }
    }

    // --- per-plane-wave K=k+G, |K|, j_l(|K|R), j_l'(|K|R), and the matching coefficients a_l,b_l ----
    std::vector<rvec3_t> K(n);
    std::vector<double>  Kmag(n);
    std::vector<std::vector<double>> a(n),b(n);
    for (size_t i=0; i<n; i++)
    {
        K[i]=B.ToCartesian(itsK+itsG[i]);
        Kmag[i]=std::sqrt(K[i]*K[i]);
        double KR=Kmag[i]*R;
        std::vector<double> jK =qchem::Math::SphericalBessel(L,KR);
        std::vector<double> jKp=qchem::Math::SphericalBesselPrime(L,KR,jK);
        a[i].assign(L+1,0.0); b[i].assign(L+1,0.0);
        for (int l=0; l<=L; l++)
        {
            double rhs1=jK[l];               // j_l(K R)
            double rhs2=Kmag[i]*jKp[l];      // K j_l'(K R)
            double W=uR[l]*udpR[l]-upR[l]*udR[l];          // Wronskian
            a[i][l]=(rhs1*udpR[l]-rhs2*udR[l])/W;
            b[i][l]=(uR[l]*rhs2-upR[l]*rhs1)/W;
        }
    }

    // --- assemble H and O = interstitial + muffin-tin -----------------------------------------------
    H=blazem::zeroH<dcmplx>(n);
    O=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            double OI;
            if (i==j) OI=1.0 - Vs/Omega;
            else { rvec3_t dG=K[i]-K[j]; double dg=std::sqrt(dG*dG);
                   OI=-(FourPi*R*R/Omega)*qchem::Math::SphericalBessel1(dg*R)/dg; }
            double KdotK=K[i]*K[j];
            double cosg=(Kmag[i]>1e-12 && Kmag[j]>1e-12) ? KdotK/(Kmag[i]*Kmag[j]) : 1.0;
            cosg=std::max(-1.0,std::min(1.0,cosg));
            std::vector<double> P=qchem::Math::LegendreP(L,cosg);
            double sphereO=0.0, sphereH=0.0;
            for (int l=0; l<=L; l++)
            {
                double ai=a[i][l], bi=b[i][l], aj=a[j][l], bj=b[j][l];
                double pref=(2*l+1)*P[l];
                sphereO += pref*( ai*aj*Suu[l] + (ai*bj+bi*aj)*Suud[l] + bi*bj*Sudud[l] );
                sphereH += pref*( ai*aj*Tuu[l] + (ai*bj+bi*aj)*Tuud[l] + bi*bj*Tudud[l] );
            }
            sphereO *= FourPi/Omega;
            sphereH *= FourPi/Omega;
            O(i,j)=dcmplx(OI + sphereO, 0.0);
            H(i,j)=dcmplx(0.5*KdotK*OI + sphereH, 0.0);
        }
}

chmat_t LAPW_IBS::MakeOverlap()     const { chmat_t H,O; Assemble(H,O); return O; }
chmat_t LAPW_IBS::MakeHamiltonian() const { chmat_t H,O; Assemble(H,O); return H; }

cvec_t LAPW_IBS::operator()(const rvec3_t& r) const
{
    size_t n=GetNumFunctions();
    double inv=1.0/sqrt(itsVolume);
    cvec_t v(n);
    for (size_t i=0; i<n; i++)
    {
        double phase=itsRecip.GetCell().ToCartesian(itsK+itsG[i])*r;
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
        rvec3_t K=itsRecip.GetCell().ToCartesian(itsK+itsG[i]);
        double phase=K*r;
        dcmplx val=dcmplx(cos(phase),sin(phase))*inv;
        g[i]=vec3_t<dcmplx>(im*K.x*val, im*K.y*val, im*K.z*val);
    }
    return g;
}

std::string LAPW_IBS::BasisSetID() const
{
    return Name()+"|nG="+std::to_string(itsG.size())+"|Rmt="+std::to_string(itsRmt)
                 +"|lmax="+std::to_string(itsLmax)+"|El="+std::to_string(itsElin);
}

std::ostream& LAPW_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " plane waves, Rmt=" << itsRmt
              << " lmax=" << itsLmax << " Elin=" << itsElin << ", " << GetSymmetry();
}

} //namespace
