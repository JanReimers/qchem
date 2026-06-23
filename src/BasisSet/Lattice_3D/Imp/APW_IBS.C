// File: BasisSet/Lattice_3D/Imp/APW_IBS.C  Augmented Plane Wave secular-matrix assembly.
module;
#include <cassert>
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.APW_IBS;
import qchem.Symmetry.Factory;   // BlochFactory
import qchem.Math;               // Pi, FourPi, Cube
import qchem.SpecialFunctions;   // SphericalBessel, SphericalBessel1, SphericalBesselPrime, LegendreP
import qchem.BasisSet.Lattice_3D.Internal.GVectors;   // BuildGs
import qchem.Blaze;              // zeroH
import qchem.Vector3D;           // dot product (operator*), norm

namespace
{
// |k+G| below this is treated as zero (the Gamma-point G=0 wave): its augmentation direction is
// irrelevant (only l=0 survives, since j_l(0)=delta_l0).
constexpr double kZeroTol = 1e-12;   // TODO: Make External
}

namespace BasisSet::Lattice_3D
{

APW_IBS::APW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
                 double Ecut, double Rmt, size_t lmax)
    : BasisSet::IrrepBasisSetImp<dcmplx>(Symmetry::BlochFactory(N,kIndex))
    , itsRecip(recip)
    , itsk(kIndex.x/static_cast<double>(N.x),
           kIndex.y/static_cast<double>(N.y),
           kIndex.z/static_cast<double>(N.z))
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
    , itsRmt(Rmt)
    , itsLmax(lmax)
{
    assert(Rmt>0.0);
    itsG = Internal::BuildGs(itsRecip, itsk, Ecut);
}

// APW has a SINGLE radial function u_l per channel (not LAPW's {u_l, udot_l} pair), so the secular
// matrix is energy-dependent -- it must be rebuilt per trial energy E (hence MakeSecular(E), not the
// constructor-time blocks of LAPW), and the sphere coupling is a scalar per l (no 2x2 quadratic form).
chmat_t APW_IBS::MakeSecular(double E) const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=GetNumFunctions();
    int    lmax=static_cast<int>(itsLmax);
    double Rmt=itsRmt, Ω=itsVolume;
    double Vsphere=FourPi*Rmt*Rmt*Rmt/3.0;
    double q=std::sqrt(2.0*E), qR=q*Rmt;

    // Per-l radial log-derivative u_l'(R)/u_l(R) = q j_l'(qR)/j_l(qR).
    rvec_t jq =SpecialFunctions::SphericalBessel(lmax,qR);
    rvec_t jqp=SpecialFunctions::SphericalBesselPrime(lmax,qR,jq);
    rvec_t logd(lmax+1,0.0);
    for (int l=0; l<=lmax; l++)
    {
        logd[l]=q*jqp[l]/jq[l];             // u_l'(R)/u_l(R)
        assert(std::isfinite(logd[l]));     // finite => not on the APW asymptote (j_l(qR)=0)
    }

    // Per-plane-wave Cartesian K=k+G, magnitude, and j_l(|K|R).
    std::vector<rvec3_t> K(n);
    rvec_t  Knorm(n);
    std::vector<rvec_t> jKR(n);
    for (size_t i=0; i<n; i++)
    {
        K[i]=B.ToCartesian(itsk+itsG[i]);
        Knorm[i]=norm(K[i]);
        jKR[i]=SpecialFunctions::SphericalBessel(lmax, Knorm[i]*Rmt);
    }

    chmat_t G=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            double KdotK=K[i]*K[j];
            // Interstitial overlap O^I = full-cell PW overlap minus the inside-sphere part.
            double overlapItstl;
            if (i==j) overlapItstl=1.0 - Vsphere/Ω;
            else { rvec3_t dG=K[i]-K[j]; double dg=norm(dG);
                   overlapItstl=-(FourPi*Rmt*Rmt/Ω)*SpecialFunctions::SphericalBessel1(dg*Rmt)/dg; }
            double interstitial=(0.5*KdotK - E)*overlapItstl;          // (1/2 K.K' - E) * O^I
            // Sphere/surface term: Sum_l (2l+1) P_l(cos g) j_l(K_i R) j_l(K_j R) u_l'(R)/u_l(R).
            double cos_γij=(Knorm[i]>kZeroTol && Knorm[j]>kZeroTol) ? KdotK/(Knorm[i]*Knorm[j]) : 1.0;
            cos_γij=std::max(-1.0,std::min(1.0,cos_γij));
            rvec_t P=SpecialFunctions::LegendreP(lmax,cos_γij);
            double sphere=0.0;
            for (int l=0; l<=lmax; l++) sphere += (2*l+1)*P[l]*jKR[i][l]*jKR[j][l]*logd[l];
            sphere *= 2.0*Pi*Rmt*Rmt/Ω;
            G(i,j)=dcmplx(interstitial+sphere, 0.0);   // single sphere at the origin => real
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
        double phase=itsRecip.GetCell().ToCartesian(itsk+itsG[i])*r; // (k+G).r  (interstitial value)
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
        rvec3_t K=itsRecip.GetCell().ToCartesian(itsk+itsG[i]);
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
