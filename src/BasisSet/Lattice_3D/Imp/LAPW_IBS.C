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
module;
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.LAPW_IBS;
import qchem.Symmetry.Factory;   // BlochFactory
import qchem.Math;               // Pi, FourPi, Cube
import qchem.SpecialFunctions;   // SphericalBessel, SphericalBessel1, SphericalBesselPrime, LegendreP
import qchem.Blaze;              // zeroH, linspace
import qchem.Vector3D;           // dot product (operator*), norm

namespace
{
// |k+G| below this is treated as zero (the Gamma-point G=0 wave): the j_l'(0) singularity is skipped
// and the augmentation direction is irrelevant (only l=0 survives, since j_l(0)=delta_l0).
constexpr double kZeroTol = 1e-12;   // TODO: Make External

// Raw radial functions on the grid: u_l=R_l, u_l', and (energy derivatives) udot_l, udot_l'.
struct RadialTable { rvec_t u,up,ud,udp; double uR,upR,udR,udpR; };

// Regular radial solution u_l=R_l and u_l' on the grid (V=-Z/r, energy E).  Z=0: analytic Bessel
// j_l(qr).  Z!=0: RK4 integration of P''=f P with P=r R_l, f(r)=l(l+1)/r^2 - 2Z/r - 2E, started from
// the r^{l+1} series at r[1] (regular at the origin).
void SolveRadial(int l,double E,double Z,const rvec_t& r,double Rmt,
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
        rvec_t jR=SpecialFunctions::SphericalBessel(l,q*Rmt);
        rvec_t jpR=SpecialFunctions::SphericalBesselPrime(l,q*Rmt,jR);
        uR=jR[l]; upR=q*jpR[l];
        return;
    }
    auto f=[&](double rr){ return l*(l+1)/(rr*rr) - 2.0*Z/rr - 2.0*E; };
    double h=Rmt/NQ;
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

// u_l, u_l', udot_l, udot_l' (udot = d/dE).  Z=0: analytic.  Z!=0: u from the ODE at E, udot by central
// finite difference of the ODE solution in E.
RadialTable BuildRadial(int l,double E,double Z,const rvec_t& r,double Rmt)
{
    int NQ=static_cast<int>(r.size())-1;
    RadialTable t;
    if (Z==0.0)
    {
        double q=std::sqrt(2.0*E);
        t.u=rvec_t(NQ+1,0.0); t.up=rvec_t(NQ+1,0.0); t.ud=rvec_t(NQ+1,0.0); t.udp=rvec_t(NQ+1,0.0);
        for (int i=1;i<=NQ;i++){ double rr=r[i],x=q*rr;
            rvec_t j=SpecialFunctions::SphericalBessel(l,x);
            rvec_t jp=SpecialFunctions::SphericalBesselPrime(l,x,j);
            double jlpp=-(2.0/x)*jp[l]-(1.0-l*(l+1)/(x*x))*j[l];
            t.u[i]=j[l]; t.up[i]=q*jp[l]; t.ud[i]=(rr/q)*jp[l]; t.udp[i]=(1.0/q)*jp[l]+rr*jlpp; }
        double xR=q*Rmt;
        rvec_t jR=SpecialFunctions::SphericalBessel(l,xR);
        rvec_t jpR=SpecialFunctions::SphericalBesselPrime(l,xR,jR);
        double jRpp=-(2.0/xR)*jpR[l]-(1.0-l*(l+1)/(xR*xR))*jR[l];
        t.uR=jR[l]; t.upR=q*jpR[l]; t.udR=(Rmt/q)*jpR[l]; t.udpR=(1.0/q)*jpR[l]+Rmt*jRpp;
        return t;
    }
    double d=1e-4;
    rvec_t u0,up0,up_,upp,um,upm; double uR0,upR0,uRp,upRp,uRm,upRm;
    SolveRadial(l,E,    Z,r,Rmt,u0,up0,uR0,upR0);
    SolveRadial(l,E+d,  Z,r,Rmt,up_,upp,uRp,upRp);
    SolveRadial(l,E-d,  Z,r,Rmt,um,upm,uRm,upRm);
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
    // Assemble the blocks once, at construction (everything below is fixed by the inputs above).
    std::vector<RadialBlock>   blocks = MuffinTinRadialBlocks();   // Act 1
    std::vector<AugmentedWave> waves  = MatchAugmentation(blocks); // Act 2
    CombineBlocks(blocks,waves);                                   // Act 3
}

// Act 1: inside the sphere, solve u_l & udot_l and integrate their overlap / kinetic / potential blocks.
std::vector<LAPW_IBS::RadialBlock> LAPW_IBS::MuffinTinRadialBlocks() const
{
    int    lmax=static_cast<int>(itsLmax);
    double Rmt=itsRmt, Elin=itsElin, Z=itsZ;
    constexpr int NQ=800;                                  // TODO: Make External (radial-quadrature points)
    double h=Rmt/NQ;
    rvec_t r=blazem::linspace(NQ+1, 0.0, Rmt);             // radial grid r[0]=0 .. r[NQ]=Rmt

    std::vector<RadialBlock> blocks(lmax+1);
    for (int l=0;l<=lmax;l++)
    {
        RadialTable t=BuildRadial(l,Elin,Z,r,Rmt);
        double Suu=0,Suud=0,Sudud=0, Kuu=0,Kuud=0,Kudud=0, Vuu=0,Vuud=0,Vudud=0;
        for (int i=1;i<=NQ;i++)                            // i=0 (r=0) contributes nothing
        {
            double w=(i==NQ ? 1.0 : (i%2 ? 4.0 : 2.0))*h/3.0;   // composite Simpson
            double rr=r[i], r2=rr*rr;
            double u=t.u[i],up=t.up[i],ud=t.ud[i],udp=t.udp[i];
            Suu   += w* u*u*r2;   Suud  += w* u*ud*r2;   Sudud += w* ud*ud*r2;
            Kuu   += w*(up*up*r2   + l*(l+1)*u*u);              // <p^2> = int(u'^2 r^2 + l(l+1)u^2)
            Kuud  += w*(up*udp*r2  + l*(l+1)*u*ud);
            Kudud += w*(udp*udp*r2 + l*(l+1)*ud*ud);
            if (Z!=0.0)                                        // <V> = int u_a (-Z/r) u_b r^2 dr
            {
                Vuu   += w*(-Z)*u*u*rr;  Vuud += w*(-Z)*u*ud*rr;  Vudud += w*(-Z)*ud*ud*rr;
            }
        }
        blocks[l].boundary =rmat2d_t(t.uR, t.udR, t.upR, t.udpR);  // [[u,udot],[u',udot']] at Rmt
        blocks[l].overlap  =rmat2d_t(Suu, Suud, Suud, Sudud);
        blocks[l].kinetic  =rmat2d_t(Kuu, Kuud, Kuud, Kudud);
        blocks[l].potential=rmat2d_t(Vuu, Vuud, Vuud, Vudud);
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

    std::vector<AugmentedWave> waves(n);
    for (size_t i=0;i<n;i++)
    {
        rvec3_t K=B.ToCartesian(itsk+itsG[i]);
        double  Knorm=norm(K);
        rvec_t  jK =SpecialFunctions::SphericalBessel(lmax,Knorm*Rmt);
        rvec_t  jKp(lmax+1,0.0);
        if (Knorm*Rmt>kZeroTol) jKp=SpecialFunctions::SphericalBesselPrime(lmax,Knorm*Rmt,jK);
        waves[i].K=K; waves[i].Knorm=Knorm; waves[i].c.resize(lmax+1);
        for (int l=0;l<=lmax;l++)
            waves[i].c[l]=boundaryInv[l]*rvec2d_t(jK[l], Knorm*jKp[l]); // K=0 -> rhs=(delta_l0,0)
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

    itsOvlp=blazem::zeroH<dcmplx>(n);
    itsKp2 =blazem::zeroH<dcmplx>(n);
    itsVnuc=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            // interstitial overlap: full-cell PW overlap (delta_ij) minus the part inside the sphere.
            double overlapItstl;
            if (i==j) overlapItstl=1.0 - Vsphere/Ω;
            else { rvec3_t dG=waves[i].K-waves[j].K; double dg=norm(dG);
                   overlapItstl=-(FourPi*Rmt*Rmt/Ω)*SpecialFunctions::SphericalBessel1(dg*Rmt)/dg; }
            double KdotK=waves[i].K*waves[j].K;
            double cos_γij=(waves[i].Knorm>kZeroTol && waves[j].Knorm>kZeroTol)
                           ? KdotK/(waves[i].Knorm*waves[j].Knorm) : 1.0;
            cos_γij=std::max(-1.0,std::min(1.0,cos_γij));
            rvec_t P=SpecialFunctions::LegendreP(lmax,cos_γij);

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
                 +"|Z="+std::to_string(itsZ);
}

std::ostream& LAPW_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " plane waves, Rmt=" << itsRmt
              << " lmax=" << itsLmax << " Elin=" << itsElin << " Z=" << itsZ << ", " << GetSymmetry();
}

} //namespace
