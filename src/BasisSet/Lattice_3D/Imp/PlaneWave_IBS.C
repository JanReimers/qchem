// File: BasisSet/Lattice_3D/Imp/PlaneWave_IBS.C  Plane-wave irrep basis set implementation.
module;
#include <algorithm>
#include <cassert>
#include <complex>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
import qchem.Symmetry.Factory;   // BlochFactory (the convenience ctor builds the Bloch irrep)
import qchem.Symmetry.BlochQN;   // Symmetry::Getk (prys k out of the abstract Bloch irrep)
import qchem.Structure;          // Atom (itsZ, itsR) + atom iteration for MakeNuclear
import qchem.Math;               // Pi, FourPi, sqrt, cos, sin, pow, Cube
import qchem.SpecialFunctions;   // LegendreP (the (2l+1)P_l angular factor)
import qchem.BasisSet.Lattice_3D.Internal.GVectors;   // BuildGs
import qchem.BasisSet.Lattice_3D.Internal.KPlusG;     // KPlusG (Cartesian k+G, |k+G|, cos gamma)
import qchem.FFT;                                      // FFT3D / NextPow2 (the XC G-space<->real transforms)
import qchem.Blaze;
import qchem.Vector3D;           // dot product (operator*) + vector arithmetic

namespace BasisSet::Lattice_3D
{

PlaneWave_IBS::PlaneWave_IBS(const ReciprocalLattice& recip, const sym_t& irrep, double Ecut)
    : BasisSet::IrrepBasisSetImp<dcmplx>(irrep)
    , itsRecip(recip)
    , itsk(Symmetry::Getk(irrep))   // the Bloch irrep IS the k-label; pry it out (mirrors Getl on atoms)
    , itsEcut(Ecut)
    // |det B| = (2 pi)^3 / |det A|, so the direct cell volume V = (2 pi)^3 / V_recip.
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
{
    itsG = Internal::BuildGs(itsRecip, itsk, Ecut);   // { G : 1/2|k+G|^2 < Ecut }
}

// Convenience: build the Bloch irrep from BZ-grid indices and delegate to the primary constructor.
PlaneWave_IBS::PlaneWave_IBS(const ReciprocalLattice& recip, const ivec3_t& N,
                             const ivec3_t& kIndex, double Ecut)
    : PlaneWave_IBS(recip, Symmetry::BlochFactory(N,kIndex), Ecut)
{}

rvec3_t PlaneWave_IBS::GetGCartesian(const ivec3_t& m) const
{
    return itsRecip.GetCell().ToCartesian(rvec3_t(m)); // B m
}

// Uniform N1xN2xN3 grid of FRACTIONAL coordinates r=(i1/N1,i2/N2,i3/N3) -- the XC real-space grid
// (uniform weight Omega/prod(N); dG.r = 2 pi dm.r_frac makes the rho(G)<->rho(r) transform a plain DFT).
std::vector<rvec3_t> PlaneWave_IBS::UniformGrid(const ivec3_t& n) const
{
    std::vector<rvec3_t> g;
    g.reserve(static_cast<size_t>(n.x)*n.y*n.z);
    for (int i1=0;i1<n.x;i1++)
        for (int i2=0;i2<n.y;i2++)
            for (int i3=0;i3<n.z;i3++)
                g.push_back(rvec3_t(i1/double(n.x), i2/double(n.y), i3/double(n.z)));
    return g;
}

namespace
{
//! Forward-transform a real field sampled on the fractional grid to its Fourier components
//! \f$\tilde V(\Delta m)=\frac1{N}\sum_r V(r)e^{-i2\pi\Delta m\cdot r}\f$ over the difference set
//! \f$\{m_i-m_j\}\f$ (the only components the matrix \f$\langle G_i|V|G_j\rangle\f$ needs).
FourierMap
ForwardDFTDiffSet(const std::vector<ivec3_t>& G, const std::vector<rvec3_t>& frac,
                  const std::vector<double>& field)
{
    FourierMap vt;
    size_t n=G.size(), Npts=frac.size();
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            ivec3_t dm=G[i]-G[j];
            if (vt.find(dm)!=vt.end()) continue;
            dcmplx s(0.0);
            for (size_t q=0;q<Npts;q++)
            {
                double ph=-2*Pi*(dm.x*frac[q].x + dm.y*frac[q].y + dm.z*frac[q].z);
                s += field[q]*dcmplx(cos(ph),sin(ph));
            }
            vt[dm]=s/double(Npts);
        }
    return vt;
}
} //anon

// Grid divisions resolving the difference set without aliasing: N > 2*(2*maxComp).
ivec3_t PlaneWave_IBS::AutoGrid() const
{
    int m=0;
    for (const ivec3_t& g : itsG)
    {
        int ax=g.x<0?-g.x:g.x, ay=g.y<0?-g.y:g.y, az=g.z<0?-g.z:g.z;
        m=std::max(m, std::max(ax, std::max(ay,az)));
    }
    int nn=4*m+1;
    return ivec3_t(nn,nn,nn);
}

// <i|f|j> weighted overlap for a real-space scalar field f: sample f on the (Cartesian) grid, forward-DFT
// to f-tilde(dm), assemble.  The XC term passes f(r)=v_xc(rho(r)); the integration is OUR business.
chmat_t PlaneWave_IBS::Overlap(const ScalarFunction<double>& f) const
{
    std::vector<rvec3_t> frac=UniformGrid(AutoGrid());
    UnitCell A=itsRecip.GetCell().MakeReciprocalCell();          // direct cell (reciprocal of the reciprocal)
    std::vector<double> field(frac.size());
    for (size_t q=0;q<frac.size();q++) field[q]=f(A.ToCartesian(frac[q]));
    auto vt=ForwardDFTDiffSet(itsG,frac,field);
    return MakePotential([&vt](const ivec3_t& dm)->dcmplx
        { auto it=vt.find(dm); return it==vt.end()?dcmplx(0.0):it->second; });
}

// Coulomb repulsion matrix + energy for a density rho: rho~(dm) via the grid DFT, then V_Coul~(dm)=4 pi
// rho~/|G|^2 (dm=0 dropped, neutralising background); E = (Omega/2) Sum_{G!=0} 4 pi |rho~|^2/|G|^2.
// From a real-space density: sample on the grid, forward-DFT to rho-tilde, then the G-space (1/r12) solve.
chmat_t PlaneWave_IBS::Repulsion(const ScalarFunction<double>& rho, double& Eh) const
{
    std::vector<rvec3_t> frac=UniformGrid(AutoGrid());
    UnitCell A=itsRecip.GetCell().MakeReciprocalCell();
    std::vector<double> field(frac.size());
    for (size_t q=0;q<frac.size();q++) field[q]=rho(A.ToCartesian(frac[q]));
    return Repulsion(ForwardDFTDiffSet(itsG,frac,field), Eh);
}

// Hartree directly from the density's G-space coefficients rho-tilde (the FFT-free Poisson solve):
//   V_H(dm) = 4 pi rho-tilde(dm)/|G|^2,   E_H = (Omega/2) Sum_{G!=0} 4 pi |rho-tilde|^2/|G|^2,  dm=0 dropped.
chmat_t PlaneWave_IBS::Repulsion(const FourierMap& rg, double& Eh) const
{
    Eh=0.0;
    for (const auto& kv : rg)
    {
        if (kv.first.x==0 && kv.first.y==0 && kv.first.z==0) continue;
        rvec3_t G=GetGCartesian(kv.first);
        Eh += FourPi*std::norm(kv.second)/(G*G);
    }
    Eh *= 0.5*itsVolume;
    return MakePotential([this,&rg](const ivec3_t& dm)->dcmplx
    {
        if (dm.x==0 && dm.y==0 && dm.z==0) return dcmplx(0.0);
        rvec3_t G=GetGCartesian(dm);
        auto it=rg.find(dm);
        return FourPi*(it==rg.end()?dcmplx(0.0):it->second)/(G*G);
    });
}

// rho-tilde(dm) = (1/Omega) Sum_{i,j: G_i-G_j=dm} D_ij.  D is Hermitian, so rho-tilde(-dm)=conj(rho-tilde(dm))
// falls out automatically (the (j,i) pair contributes conj(D_ij) at -dm).  One O(n^2) pass, no grid.
FourierMap PlaneWave_IBS::MakeFourierDensity(const chmat_t& D) const
{
    FourierMap rg;
    size_t n=GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
            rg[itsG[i]-itsG[j]] += D(i,j);
    for (auto& kv : rg) kv.second /= itsVolume;
    return rg;
}

// AutoGrid divisions rounded up to powers of two -- radix-2 FFT for the XC route.  A larger grid still
// resolves the difference set without aliasing, so it is at worst slightly more accurate.
ivec3_t PlaneWave_IBS::FFTGrid() const
{
    ivec3_t a=AutoGrid();
    return ivec3_t(int(qchem::FFT::NextPow2(a.x)), int(qchem::FFT::NextPow2(a.y)), int(qchem::FFT::NextPow2(a.z)));
}

// rho(r) on the FFT grid = inverse FFT of rho-tilde: rho(r_j) = Sum_dm rho-tilde(dm) e^{+i2pi dm.j/N}.
// rho-tilde is the physical coefficient (already /Omega), so the inverse FFT takes NO 1/N normalization.
rvec_t PlaneWave_IBS::RhoOnGrid(const FourierMap& rho) const
{
    ivec3_t N=FFTGrid();
    size_t Npts=size_t(N.x)*N.y*N.z;
    cvec_t g(Npts, dcmplx(0.0));
    for (const auto& kv : rho)
    {
        int i0=((kv.first.x%N.x)+N.x)%N.x, i1=((kv.first.y%N.y)+N.y)%N.y, i2=((kv.first.z%N.z)+N.z)%N.z;
        g[(size_t(i0)*N.y+i1)*N.z+i2]=kv.second;
    }
    cvec_t rr=qchem::FFT::FFT3D(g, N, +1);
    rvec_t out(Npts);
    for (size_t i=0;i<Npts;i++) out[i]=std::real(dcmplx(rr[i]));
    return out;
}

// Forward-FFT a real-space grid field to its G-space coefficients Vtilde(dm)=(1/Npts) FFT[V], stored as a
// FourierMap over the basis difference set {m_i-m_j, j>=i} -- the keys the assembly will query.  The
// potential analogue of MakeFourierDensity (which produces rho-tilde from the density matrix).
FourierMap PlaneWave_IBS::ForwardGrid(const rvec_t& V) const
{
    ivec3_t N=FFTGrid();
    size_t Npts=size_t(N.x)*N.y*N.z;
    assert(V.size()==Npts);
    cvec_t g(Npts, dcmplx(0.0));
    for (size_t i=0;i<Npts;i++) g[i]=dcmplx(V[i]);
    cvec_t Vt=qchem::FFT::FFT3D(g, N, -1);
    FourierMap out;
    size_t n=GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            ivec3_t dm=itsG[i]-itsG[j];
            if (out.find(dm)==out.end())
            {
                int i0=((dm.x%N.x)+N.x)%N.x, i1=((dm.y%N.y)+N.y)%N.y, i2=((dm.z%N.z)+N.z)%N.z;
                out[dm]=Vt[(size_t(i0)*N.y+i1)*N.z+i2]/double(Npts);
            }
        }
    return out;
}

// <i|V|j> = Vtilde(m_i-m_j): assemble directly from the G-space coefficients (no kernel -- the overlap
// 3-centre is the delta).  The XC sibling of Repulsion (which folds in 4pi/G^2).
chmat_t PlaneWave_IBS::Overlap(const FourierMap& Vt) const
{
    return MakePotential([&](const ivec3_t& dm)->dcmplx
        { auto it=Vt.find(dm); return it==Vt.end()?dcmplx(0.0):it->second; });
}

// Real-space grid -> matrix: forward-FFT then assemble.
chmat_t PlaneWave_IBS::Overlap(const rvec_t& V) const
{
    return Overlap(ForwardGrid(V));
}

// integral f d3r on the FFT grid: uniform quadrature, weight Omega/Npts.
double PlaneWave_IBS::Integral(const rvec_t& f) const
{
    return blazem::sum(f)*itsVolume/double(f.size());
}

// Scalar integral integral f d3r over the cell: uniform-grid quadrature (weight Omega/Npts).
double PlaneWave_IBS::Integral(const ScalarFunction<double>& f) const
{
    std::vector<rvec3_t> frac=UniformGrid(AutoGrid());
    UnitCell A=itsRecip.GetCell().MakeReciprocalCell();
    double s=0.0;
    for (const rvec3_t& p : frac) s += f(A.ToCartesian(p));
    return s*itsVolume/double(frac.size());
}

// Energy of the DROPPED G=0 local-potential component: E_alpha = (N/Omega) Sum_a alpha_a, with
// alpha_a = the local model's finite G->0 limit (FormFactorG0 = integral[V_loc+Z/r]).  itsVolume = Omega.
// Only the LOCAL part has a G=0 alignment (the separable nonlocal is short-ranged); a bare-Coulomb model
// has alpha=0.  This is an ENERGY-only term -- the G=0 potential never entered the matrix
// (MakeLocalPotential drops dG=0), so it is added to the total energy by the external term, which owns
// the model and supplies it here (Omega and the cell geometry stay the basis's business).
double PlaneWave_IBS::ExternalG0Energy(const Structure* cl, const LocalPotential& loc, double numElectrons) const
{
    double sumAlpha=0.0;
    for (Atom* a : *cl) sumAlpha += loc.FormFactorG0(a->itsZ);
    return (numElectrons/itsVolume)*sumAlpha;
}

// Plane waves are orthonormal over the cell: <G|G'> = delta_{GG'}.
chmat_t PlaneWave_IBS::MakeOverlap() const
{
    size_t n=GetNumFunctions();
    chmat_t S=blazem::zeroH<dcmplx>(n);   // hmat_t(n) does NOT zero the off-diagonals
    for (size_t i=0; i<n; i++) S(i,i)=1.0;
    return S;
}

// <p^2> = <-nabla^2> building block (NO 1/2 -- the Hamiltonian applies it).  For a plane wave
// -nabla^2 e^{i(k+G).r} = |k+G|^2 e^{i(k+G).r}, so the matrix is diagonal in |k+G|^2.
chmat_t PlaneWave_IBS::MakeKinetic() const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=GetNumFunctions();
    chmat_t S=blazem::zeroH<dcmplx>(n);   // off-diagonals are exactly zero
    for (size_t i=0; i<n; i++)
    {
        double kG=B.GetDistance(itsk+itsG[i]); // |k+G|
        S(i,i)=kG*kG;
    }
    return S;
}

// Bare-Coulomb electron-nucleus attraction in reciprocal space:
//   <G|V|G'> = V(dG) = -(4 pi / Omega) Sum_a Z_a e^{-i dG.tau_a} / |dG|^2,   dG = G-G' != 0.
// The dG=0 term (the divergent G=0 Coulomb component) is dropped -- the conventional uniform
// neutralising background; it contributes only a finite per-cell shift that -> 0 as the cell grows.
// Bare nuclear Coulomb is just the local potential with the BareCoulomb form factor.
chmat_t PlaneWave_IBS::MakeNuclear(const Structure* cl) const
{
    return MakeLocalPotential(cl, BareCoulomb());
}

// V(dG) = (1/Omega) Sum_a v(Z_a,|dG|^2) e^{-i dG.tau_a}, dG=0 dropped (neutralising background).
// The result is Hermitian: V(-dG) = conj(V(dG)) since the structure factor conjugates under dG -> -dG
// (the real form factor v is even).  Filling the upper triangle of a HermitianMatrix auto-sets the
// lower as the conjugate, so off-origin / multi-atom cells (complex phases) are handled correctly.
chmat_t PlaneWave_IBS::MakeLocalPotential(const Structure* cl, const LocalPotential& v) const
{
    const UnitCell& B=itsRecip.GetCell();
    return MakePotential([&](const ivec3_t& dm)->dcmplx
    {
        if (dm.x==0 && dm.y==0 && dm.z==0) return dcmplx(0.0); // drop dG=0
        rvec3_t dG=B.ToCartesian(rvec3_t(dm));                  // dG = B.dm (Cartesian)
        double g2=dG*dG;
        dcmplx acc(0.0);                                        // (form factor) x (structure factor)
        for (Atom* a : *cl) acc += v.FormFactor(a->itsZ,g2)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
        return acc/itsVolume;
    });
}

// V_NL(G,G') = (1/Omega) Sum_a e^{-i(G-G').tau_a} Sum_p (2l_p+1) P_{l_p}(cos gamma) betã_p(|k+G|) D_p
//              betã_p(|k+G'|),  gamma = angle(k+G, k+G').
// The angular factor (2l+1)P_l(cos gamma) is the addition-theorem sum Sum_m Y_lm(q^)Y*_lm(q'^) over the
// projector's m -- the same structure the APW/LAPW sphere terms use; l=0 gives P_0=1 (the s-channel).
// (The k-point phases of <k+G|beta>, <beta|k+G'> cancel, leaving the structure-factor phase.)
// Per atom & projector & m this is rank-1: |beta> D <beta|.  Hermitian; real for atoms at the origin.
chmat_t PlaneWave_IBS::MakeSeparablePotential(const Structure* cl, const SeparablePotential& v) const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=GetNumFunctions();
    Internal::KPlusG kg(B, itsk, itsG);             // Cartesian k+G, |k+G|, and cos(gamma)

    int maxL=0;                                     // highest projector channel present
    for (Atom* a : *cl)
        for (size_t p=0; p<v.NumProjectors(a->itsZ); p++)
            maxL=std::max(maxL, v.AngularMomentum(a->itsZ,p));

    chmat_t V=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            rvec3_t dG=B.ToCartesian(rvec3_t(itsG[i]-itsG[j]));
            rvec_t P=SpecialFunctions::LegendreP(maxL, kg.CosGamma(i,j));
            dcmplx acc(0.0);
            for (Atom* a : *cl)
            {
                double s=0.0;                       // Sum_p (2l+1)P_l(cos gamma) betã_p(q_i) D_p betã_p(q_j)
                for (size_t p=0; p<v.NumProjectors(a->itsZ); p++)
                {
                    int l=v.AngularMomentum(a->itsZ,p);
                    s += (2*l+1)*P[l] * v.Projector(a->itsZ,p,kg.Norm(i))*v.Coefficient(a->itsZ,p)
                                       *v.Projector(a->itsZ,p,kg.Norm(j));
                }
                acc += s*std::exp(dcmplx(0.0,-(dG*a->itsR)));
            }
            V(i,j)=acc/itsVolume;
        }
    return V;
}

// <G|V|G'> = Vtilde(m(G) - m(G')).  Fill the upper triangle; HermitianMatrix mirrors the conjugate.
chmat_t PlaneWave_IBS::MakePotential(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    size_t n=GetNumFunctions();
    chmat_t V=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
            V(i,j)=Vtilde(itsG[i]-itsG[j]);
    return V;
}

cvec_t PlaneWave_IBS::operator()(const rvec3_t& r) const
{
    size_t n=GetNumFunctions();
    double invSqrtV=1.0/sqrt(itsVolume);
    cvec_t v(n);
    for (size_t i=0; i<n; i++)
    {
        double phase=GetGCartesian(itsG[i])*r + itsRecip.GetCell().ToCartesian(itsk)*r; // (k+G).r
        v[i]=dcmplx(cos(phase),sin(phase))*invSqrtV;
    }
    return v;
}

cvec3vec_t PlaneWave_IBS::Gradient(const rvec3_t& r) const
{
    // grad e^{i(k+G).r}/sqrt(V) = i(k+G) e^{i(k+G).r}/sqrt(V).
    const dcmplx im(0.0,1.0);
    size_t n=GetNumFunctions();
    double invSqrtV=1.0/sqrt(itsVolume);
    rvec3_t kCart=itsRecip.GetCell().ToCartesian(itsk);
    cvec3vec_t g(n);
    for (size_t i=0; i<n; i++)
    {
        rvec3_t kG=kCart+GetGCartesian(itsG[i]); // k+G (Cartesian)
        double phase=kG*r;
        dcmplx val=dcmplx(cos(phase),sin(phase))*invSqrtV;
        g[i]=vec3_t<dcmplx>(im*kG.x*val, im*kG.y*val, im*kG.z*val);
    }
    return g;
}

std::string PlaneWave_IBS::BasisSetID() const
{
    return Name()+"|k="+std::to_string(itsk.x)+","+std::to_string(itsk.y)+","+std::to_string(itsk.z)
                 +"|Ecut="+std::to_string(itsEcut)
                 +"|nG="+std::to_string(itsG.size());
}

std::ostream& PlaneWave_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " plane waves, "
              << GetSymmetry();
}

} //namespace
