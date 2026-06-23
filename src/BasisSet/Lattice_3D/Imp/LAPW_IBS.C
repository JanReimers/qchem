// File: BasisSet/Lattice_3D/Imp/LAPW_IBS.C  LAPW Hamiltonian/overlap/nuclear assembly.
//
// The constructor assembles the three Orbital_1E_IBS blocks (overlap, kinetic <p^2>, nuclear <V>) in
// three physics acts:
//   Act 1  MuffinTinRadialBlocks  -- solve the radial functions u_l, udot_l inside the sphere and form
//                                     their 2x2 overlap / kinetic / potential blocks in the {u,udot} basis;
//   Act 2  MatchAugmentation      -- for each plane wave, match value AND slope at the sphere boundary
//                                     (a 2x2 solve) to get the augmentation coefficients c_l=(a_l,b_l);
//   Act 3  CombineBlocks          -- H_ij = interstitial (plane wave) + muffin-tin, where the muffin-tin
//                                     part is sum_l (2l+1) P_l(cos g) * c_i . (block_l . c_j).
//
// The radial state at a grid point is carried as a single 2x2:  Phi = [[u, udot],[u', udot']]  -- rows are
// value and slope, columns are u and udot.  Phi at r=Rmt is the boundary matrix that solves the matching;
// the radial integral blocks are Simpson integrals of outer products of its value/slope rows.
module;
#include <cassert>
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.LAPW_IBS;
import qchem.Symmetry.Factory;   // BlochFactory
import qchem.Math;               // Pi, FourPi, Cube
import qchem.SpecialFunctions;   // SphericalBessel, SphericalBessel1, SphericalBesselPrime, LegendreP
import qchem.BasisSet.Lattice_3D.Internal.GVectors;   // BuildGs
import qchem.BasisSet.Lattice_3D.Internal.KPlusG;     // KPlusG (Cartesian k+G, |k+G|, cos gamma) + kZeroTol
import qchem.Blaze;              // zeroH
import qchem.Vector3D;           // dot product (operator*), norm

namespace
{
using BasisSet::Lattice_3D::Internal::kZeroTol;  // |k+G| below this -> skip the j_l'(0) singularity

// integral f(r) dr over a LOGARITHMIC grid (r uniform in x=ln r): composite Simpson in x with the
// dr = r dx Jacobian, so integral f dr = sum_i w_i^x f(r_i) r_i.  The log grid clusters points near the
// origin, which keeps the centrifugal stiffness h^2*(l(l+1)/r^2) ~ dx^2*l(l+1) bounded for all l (the
// uniform grid was stiff for l>=1).  Generic in f's value type (here a 2x2 rmat2d_t); f(i)=integrand at r[i].
template <class F> auto RadialIntegral(const rvec_t& r, F f)
{
    int    NQ=static_cast<int>(r.size())-1;
    double dx=std::log(r[1]/r[0]);                           // uniform step in x=ln r
    auto acc=(dx/3.0)*r[0]*f(0);                             // i=0 endpoint (+ dr/dx = r Jacobian)
    for (int i=1;i<=NQ;i++)
    {
        double w=(i==NQ ? 1.0 : (i%2 ? 4.0 : 2.0))*dx/3.0;
        acc += w*r[i]*f(i);
    }
    return acc;
}

// The radial state on the grid: Phi[i] = [[u, udot],[u', udot']] at r[i] (rows value/slope, cols u/udot).
// boundary = Phi at r=Rmt (== Phi.back()), the matrix that solves the value+slope matching.
struct RadialTable { std::vector<rmat2d_t> Phi; rmat2d_t boundary; };

// (value, slope) = (u_l(r_i), u_l'(r_i)) of the regular radial solution of -1/2 grad^2 - Znuc/r at energy E,
// on the logarithmic grid r (r[0]=rmin>0, r[NQ]=Rmt).  Znuc=0: analytic Bessel j_l(qr).  Znuc!=0:
// VARIABLE-step RK4 of P''=f P (P=r R_l), f=l(l+1)/r^2 - 2 Znuc/r - 2E, from the r^{l+1} series at r[0];
// the log grid's shrinking step near the origin keeps the centrifugal term stable for l>=1.
std::vector<rvec2d_t> SolveRadial(int l,double E,double Znuc,const rvec_t& r)
{
    int NQ=static_cast<int>(r.size())-1;
    std::vector<rvec2d_t> us(NQ+1, rvec2d_t(0,0));
    if (Znuc==0.0)
    {
        assert(E>0.0);                                       // free radial solution j_l(sqrt(2E) r)
        double q=std::sqrt(2.0*E);
        for (int i=0;i<=NQ;i++){ double x=q*r[i];
            rvec_t j =SpecialFunctions::SphericalBessel(l,x);
            rvec_t jp=SpecialFunctions::SphericalBesselPrime(l,x,j);
            us[i]=rvec2d_t(j[l], q*jp[l]); }
        return us;
    }
    auto f=[&](double rr){ return l*(l+1)/(rr*rr) - 2.0*Znuc/rr - 2.0*E; };
    double P=std::pow(r[0],l+1), Q=(l+1)*std::pow(r[0],l);    // u ~ r^l series start at r[0]=rmin
    us[0]=rvec2d_t(P/r[0], Q/r[0]-P/(r[0]*r[0]));
    for (int i=0;i<NQ;i++)
    {
        double rr=r[i], h=r[i+1]-r[i];                       // variable step on the log grid
        double k1P=Q,           k1Q=f(rr)*P;
        double k2P=Q+0.5*h*k1Q, k2Q=f(rr+0.5*h)*(P+0.5*h*k1P);
        double k3P=Q+0.5*h*k2Q, k3Q=f(rr+0.5*h)*(P+0.5*h*k2P);
        double k4P=Q+h*k3Q,     k4Q=f(rr+h)*(P+h*k3P);
        P+=h/6.0*(k1P+2*k2P+2*k3P+k4P);
        Q+=h/6.0*(k1Q+2*k2Q+2*k3Q+k4Q);
        double rn=r[i+1];
        us[i+1]=rvec2d_t(P/rn, Q/rn-P/(rn*rn));
    }
    return us;
}

// Build Phi[i] = [[u, udot],[u', udot']] on the log grid.  Znuc=0: analytic.  Znuc!=0: u from the ODE
// at E, udot by central finite difference of the ODE solution in E.
RadialTable BuildRadial(int l,double E,double Znuc,const rvec_t& r)
{
    int NQ=static_cast<int>(r.size())-1;
    std::vector<rmat2d_t> Phi(NQ+1);
    if (Znuc==0.0)
    {
        assert(E>0.0);                                       // free radial solution j_l(sqrt(2E) r)
        double q=std::sqrt(2.0*E);
        for (int i=0;i<=NQ;i++){ double rr=r[i],x=q*rr;
            rvec_t j =SpecialFunctions::SphericalBessel(l,x);
            rvec_t jp=SpecialFunctions::SphericalBesselPrime(l,x,j);
            double jlpp=-(2.0/x)*jp[l]-(1.0-l*(l+1)/(x*x))*j[l];
            double u=j[l], up=q*jp[l], ud=(rr/q)*jp[l], udp=(1.0/q)*jp[l]+rr*jlpp;
            Phi[i]=rmat2d_t(u, ud, up, udp); }                 // [[u,udot],[u',udot']]
    }
    else
    {
        double d=1e-4;
        std::vector<rvec2d_t> uE=SolveRadial(l,E,  Znuc,r);
        std::vector<rvec2d_t> uP=SolveRadial(l,E+d,Znuc,r);
        std::vector<rvec2d_t> uM=SolveRadial(l,E-d,Znuc,r);
        for (int i=0;i<=NQ;i++)
        {
            rvec2d_t udot=(uP[i]-uM[i])/(2*d);                 // (udot, udot')
            Phi[i]=rmat2d_t(uE[i].x, udot.x, uE[i].y, udot.y); // [[u,udot],[u',udot']]
        }
    }
    rmat2d_t boundary=Phi[NQ];                             // grab before moving Phi out
    return RadialTable{ std::move(Phi), boundary };
}
} // anon namespace

namespace BasisSet::Lattice_3D
{

LAPW_IBS::LAPW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
                   double Ecut, double Rmt, size_t lmax, double Elin, double Znuc)
    : BasisSet::IrrepBasisSetImp<dcmplx>(Symmetry::BlochFactory(N,kIndex))
    , itsRecip(recip)
    , itsk(kIndex.x/static_cast<double>(N.x),
           kIndex.y/static_cast<double>(N.y),
           kIndex.z/static_cast<double>(N.z))
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
    , itsRmt(Rmt)
    , itsLmax(lmax)
    , itsElin(Elin)
    , itsZnuc(Znuc)
{
    assert(Rmt>0.0);
    itsG = Internal::BuildGs(itsRecip, itsk, Ecut);
    // Assemble the blocks once, at construction (everything below is fixed by the inputs above).
    std::vector<RadialBlock>   blocks = MuffinTinRadialBlocks();   // Act 1
    std::vector<AugmentedWave> waves  = MatchAugmentation(blocks); // Act 2
    CombineBlocks(blocks,waves);                                   // Act 3
}

// Act 1: inside the sphere, solve u_l & udot_l and integrate their overlap / kinetic / potential blocks.
// Each block is a Simpson integral of an outer product of Phi's value row v=(u,udot) and/or slope row
// s=(u',udot'):  overlap = int r^2 v(x)v;  kinetic = int [r^2 s(x)s + l(l+1) v(x)v];  V = int (-Znuc r) v(x)v.
std::vector<LAPW_IBS::RadialBlock> LAPW_IBS::MuffinTinRadialBlocks() const
{
    int    lmax=static_cast<int>(itsLmax);
    double Rmt=itsRmt, Elin=itsElin, Znuc=itsZnuc;
    constexpr int    NQ=800;                               // TODO: Make External (radial-quadrature points)
    constexpr double rmin=1e-5;                            // TODO: Make External (log-grid inner radius, Bohr)
    rvec_t r(NQ+1);                                        // logarithmic grid r[0]=rmin .. r[NQ]=Rmt
    double dx=std::log(Rmt/rmin)/NQ;
    for (int i=0;i<=NQ;i++) r[i]=rmin*std::exp(i*dx);

    std::vector<RadialBlock> blocks(lmax+1);
    for (int l=0;l<=lmax;l++)
    {
        RadialTable t=BuildRadial(l,Elin,Znuc,r);
        rmat2d_t overlap  =RadialIntegral(r,[&](int i){ rvec2d_t v=t.Phi[i].GetRow(1);
                                                        return r[i]*r[i]*Outer(v,v); });
        rmat2d_t kinetic  =RadialIntegral(r,[&](int i){ rvec2d_t v=t.Phi[i].GetRow(1), s=t.Phi[i].GetRow(2);
                                                        return r[i]*r[i]*Outer(s,s) + double(l*(l+1))*Outer(v,v); });
        rmat2d_t potential=RadialIntegral(r,[&](int i){ rvec2d_t v=t.Phi[i].GetRow(1);
                                                        return (-Znuc*r[i])*Outer(v,v); });
        blocks[l]=RadialBlock{ t.boundary, overlap, kinetic, potential };
    }
    return blocks;
}

// Act 2: augment each plane wave -- match value AND slope at Rmt.  The two matching conditions are the
// 2x2 system  boundary . (a_l,b_l) = (j_l(KR), K j_l'(KR)),  so c_l = boundary^{-1} . rhs.
std::vector<LAPW_IBS::AugmentedWave> LAPW_IBS::MatchAugmentation(const std::vector<RadialBlock>& blocks) const
{
    const UnitCell& B=itsRecip.GetCell();
    int    lmax=static_cast<int>(itsLmax);
    double Rmt=itsRmt;
    size_t n=GetNumFunctions();

    std::vector<rmat2d_t> boundaryInv(lmax+1);             // depends only on l, so invert once
    for (int l=0;l<=lmax;l++) boundaryInv[l]=Invert(blocks[l].boundary);

    Internal::KPlusG kg(B, itsk, itsG);                    // Cartesian k+G, |k+G|
    std::vector<AugmentedWave> waves;
    waves.reserve(n);
    for (size_t i=0;i<n;i++)
    {
        double  Knorm=kg.Norm(i);
        rvec_t  jK =SpecialFunctions::SphericalBessel(lmax,Knorm*Rmt);
        rvec_t  jKp(lmax+1,0.0);
        if (Knorm*Rmt>kZeroTol) jKp=SpecialFunctions::SphericalBesselPrime(lmax,Knorm*Rmt,jK);
        std::vector<rvec2d_t> c(lmax+1);
        for (int l=0;l<=lmax;l++)
        {
            c[l]=boundaryInv[l]*rvec2d_t(jK[l], Knorm*jKp[l]); // K=0 -> rhs=(delta_l0,0)
            assert(std::isfinite(c[l].x) && std::isfinite(c[l].y)); // matching well-posed (no NaN/inf)
        }
        waves.push_back(AugmentedWave(std::move(c)));
    }
    return waves;
}

// Act 3: H_ij = interstitial (plane wave) + muffin-tin.  The muffin-tin part of each operator is an
// angular sum of the matching-coefficient quadratic form  c_i . (block_l . c_j)  of its radial block.
void LAPW_IBS::CombineBlocks(const std::vector<RadialBlock>& blocks, const std::vector<AugmentedWave>& waves)
{
    int    lmax=static_cast<int>(itsLmax);
    double Rmt=itsRmt, Ω=itsVolume;
    double Vsphere=FourPi*Rmt*Rmt*Rmt/3.0;
    size_t n=GetNumFunctions();

    Internal::KPlusG kg(itsRecip.GetCell(), itsk, itsG);   // Cartesian k+G, |k+G|, cos gamma
    itsOvlp=blazem::zeroH<dcmplx>(n);
    itsKp2 =blazem::zeroH<dcmplx>(n);
    itsVnuc=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            // interstitial overlap: full-cell PW overlap (delta_ij) minus the part inside the sphere.
            double overlapItstl;
            if (i==j) overlapItstl=1.0 - Vsphere/Ω;
            else { rvec3_t dG=kg.K(i)-kg.K(j); double dg=norm(dG);
                   overlapItstl=-(FourPi*Rmt*Rmt/Ω)*SpecialFunctions::SphericalBessel1(dg*Rmt)/dg; }
            double KdotK=kg.K(i)*kg.K(j);
            rvec_t P=SpecialFunctions::LegendreP(lmax, kg.CosGamma(i,j));

            double sphereO=0.0, sphereK=0.0, sphereV=0.0;
            for (int l=0;l<=lmax;l++)
            {
                rvec2d_t ci=waves[i].c[l], cj=waves[j].c[l];
                double w=(2*l+1)*P[l];                          // angular weight (addition theorem)
                sphereO += w * (ci*(blocks[l].overlap  *cj));   // c_i . (block . c_j)
                sphereK += w * (ci*(blocks[l].kinetic  *cj));
                sphereV += w * (ci*(blocks[l].potential*cj));
            }
            itsOvlp(i,j)=dcmplx(overlapItstl + sphereO*FourPi/Ω, 0.0);
            itsKp2 (i,j)=dcmplx(KdotK*overlapItstl + sphereK*FourPi/Ω, 0.0); // <p^2> (NO 1/2)
            itsVnuc(i,j)=dcmplx(sphereV*FourPi/Ω, 0.0);                      // interstitial V=0
        }
}

chmat_t LAPW_IBS::MakeOverlap() const { return itsOvlp; }
chmat_t LAPW_IBS::MakeKinetic() const { return itsKp2; }
chmat_t LAPW_IBS::MakeNuclear(const Structure*) const { return itsVnuc; }

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
                 +"|Z="+std::to_string(itsZnuc);
}

std::ostream& LAPW_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " plane waves, Rmt=" << itsRmt
              << " lmax=" << itsLmax << " Elin=" << itsElin << " Z=" << itsZnuc << ", " << GetSymmetry();
}

} //namespace
