// File: BasisSet/Lattice_3D/Imp/APW_IBS.C  Augmented Plane Wave secular-matrix assembly.
module;
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.APW_IBS;
import qchem.Symmetry.Factory;   // BlochFactory
import qchem.Math;               // Pi, FourPi, Cube + qchem::Math special functions (Bessel, Legendre)
import qchem.Blaze;              // zeroH
import qchem.Vector3D;           // dot product (operator*)

namespace BasisSet::Lattice_3D
{

APW_IBS::APW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
                 double Ecut, double Rmt, size_t lmax)
    : BasisSet::IrrepBasisSetImp<dcmplx>(Symmetry::BlochFactory(N,kIndex))
    , itsRecip(recip)
    , itsK(kIndex.x/static_cast<double>(N.x),
           kIndex.y/static_cast<double>(N.y),
           kIndex.z/static_cast<double>(N.z))
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
    , itsRmt(Rmt)
    , itsLmax(lmax)
{
    const UnitCell& B=itsRecip.GetCell();
    double Gmax=sqrt(2*Ecut)+B.GetDistance(itsK);
    for (const ivec3_t& m : itsRecip.GetGVectors(Gmax))
    {
        double kG=B.GetDistance(itsK+m);
        if (0.5*kG*kG < Ecut) itsG.push_back(m);
    }
}

chmat_t APW_IBS::MakeSecular(double E) const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=GetNumFunctions();
    int    L=static_cast<int>(itsLmax);
    double R=itsRmt, Omega=itsVolume;
    double Vs=FourPi*R*R*R/3.0;
    double q=std::sqrt(2.0*E), qR=q*R;

    // Per-l radial log-derivative u_l'(R)/u_l(R) = q j_l'(qR)/j_l(qR).
    std::vector<double> jq =qchem::Math::SphericalBessel(L,qR);
    std::vector<double> jqp=qchem::Math::SphericalBesselPrime(L,qR,jq);
    std::vector<double> logd(L+1,0.0);
    for (int l=0; l<=L; l++) logd[l]=q*jqp[l]/jq[l];

    // Per-plane-wave Cartesian K=k+G, magnitude, and j_l(|K|R).
    std::vector<rvec3_t> K(n);
    std::vector<double>  Kmag(n);
    std::vector<std::vector<double>> jKR(n);
    for (size_t i=0; i<n; i++)
    {
        K[i]=B.ToCartesian(itsK+itsG[i]);
        Kmag[i]=std::sqrt(K[i]*K[i]);
        jKR[i]=qchem::Math::SphericalBessel(L, Kmag[i]*R);
    }

    chmat_t G=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            double KdotK=K[i]*K[j];
            // Interstitial overlap O^I = full-cell PW overlap minus the inside-sphere part.
            double OI;
            if (i==j) OI=1.0 - Vs/Omega;
            else { rvec3_t dG=K[i]-K[j]; double dg=std::sqrt(dG*dG); OI=-(FourPi*R*R/Omega)*qchem::Math::SphericalBessel1(dg*R)/dg; }
            double gamma=(0.5*KdotK - E)*OI;
            // Sphere/surface term: Sum_l (2l+1) P_l(cos g) j_l(K_i R) j_l(K_j R) u_l'(R)/u_l(R).
            double cosg=(Kmag[i]>1e-12 && Kmag[j]>1e-12) ? KdotK/(Kmag[i]*Kmag[j]) : 1.0;
            cosg=std::max(-1.0,std::min(1.0,cosg));
            std::vector<double> P=qchem::Math::LegendreP(L,cosg);
            double sphere=0.0;
            for (int l=0; l<=L; l++) sphere += (2*l+1)*P[l]*jKR[i][l]*jKR[j][l]*logd[l];
            sphere *= 2.0*Pi*R*R/Omega;
            G(i,j)=dcmplx(gamma+sphere, 0.0);   // single sphere at the origin => real
        }
    return G;
}

cvec_t APW_IBS::operator()(const rvec3_t& r) const
{
    size_t n=GetNumFunctions();
    double inv=1.0/sqrt(itsVolume);
    cvec_t v(n);
    for (size_t i=0; i<n; i++)
    {
        double phase=itsRecip.GetCell().ToCartesian(itsK+itsG[i])*r; // (k+G).r  (interstitial value)
        v[i]=dcmplx(cos(phase),sin(phase))*inv;
    }
    return v;
}

cvec3vec_t APW_IBS::Gradient(const rvec3_t& r) const
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

std::string APW_IBS::BasisSetID() const
{
    return Name()+"|nG="+std::to_string(itsG.size())+"|Rmt="+std::to_string(itsRmt)
                 +"|lmax="+std::to_string(itsLmax);
}

std::ostream& APW_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " plane waves, Rmt=" << itsRmt
              << " lmax=" << itsLmax << ", " << GetSymmetry();
}

} //namespace
