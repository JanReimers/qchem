// File: BasisSet/Lattice_3D/Imp/PlaneWave_IBS.C  Plane-wave irrep basis set implementation.
//
// Grid-geometry methods (op(r), overlap/kinetic, the {G} set, MakePotential, the FFT grid) live in the
// shared PW_Evaluator (this basis IS-A one, reached through the EPW_* mixins).  What remains here is the
// orbital-only, atom/model-driven assembly: the density-driven G-space Hartree/XC route and the external
// pseudopotential.  Those read the shared grid data through the evaluator accessors (Gs(), Volume(),
// Recip(), kFrac(), GetGCartesian, MakePotential, the FFT-grid helpers).
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
import qchem.BasisSet.Lattice_3D.PlaneWaveFit_IBS;   // the auxiliary density-fit basis CreateCDFitBasisSet builds
import qchem.Symmetry.Factory;   // BlochFactory (the convenience ctor builds the Bloch irrep)
import qchem.Symmetry.Lattice_3D.BlochQN;   // Symmetry::Lattice_3D::Getk (prys k out of the abstract Bloch irrep)
import qchem.Structure;          // Atom (itsZ, itsR) + atom iteration for MakeNuclear
import qchem.Math;               // Pi, FourPi, cos, sin
import qchem.SpecialFunctions;   // LegendreP (the (2l+1)P_l angular factor)
import qchem.BasisSet.Lattice_3D.Internal.KPlusG;     // KPlusG (Cartesian k+G, |k+G|, cos gamma)
import qchem.FFT;                                      // FFT3D (the XC G-space<->real transforms)
import qchem.Blaze;
import qchem.Vector3D;           // dot product (operator*) + vector arithmetic

namespace qchem::BasisSet::Lattice_3D
{

PlaneWave_IBS::PlaneWave_IBS(const ReciprocalLattice& recip, const sym_t& irrep, double Ecut)
    : BasisSet::IrrepBasisSetImp<dcmplx>(irrep)
    , PW_Evaluator(recip, Symmetry::Lattice_3D::Getk(irrep), Ecut) // the Bloch irrep IS the k-label; pry it out
{}

// Convenience: build the Bloch irrep from BZ-grid indices and delegate to the primary constructor.
PlaneWave_IBS::PlaneWave_IBS(const ReciprocalLattice& recip, const ivec3_t& N,
                             const ivec3_t& kIndex, double Ecut)
    : PlaneWave_IBS(recip, Symmetry::BlochFactory(N,kIndex), Ecut)
{}

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

// <i|f|j> weighted overlap for a real-space scalar field f: sample f on the (Cartesian) grid, forward-DFT
// to f-tilde(dm), assemble.  The XC term passes f(r)=v_xc(rho(r)); the integration is OUR business.
chmat_t PlaneWave_IBS::Overlap(const ScalarFunction<double>& f) const
{
    std::vector<rvec3_t> frac=UniformGrid(AutoGrid());
    UnitCell A=Recip().GetCell().MakeReciprocalCell();          // direct cell (reciprocal of the reciprocal)
    std::vector<double> field(frac.size());
    for (size_t q=0;q<frac.size();q++) field[q]=f(A.ToCartesian(frac[q]));
    auto vt=ForwardDFTDiffSet(Gs(),frac,field);
    return MakePotential([&vt](const ivec3_t& dm)->dcmplx
        { auto it=vt.find(dm); return it==vt.end()?dcmplx(0.0):it->second; });
}

// Coulomb repulsion matrix + energy for a density rho: sample on the grid, forward-DFT to rho-tilde, then
// the G-space (1/r12) solve.
chmat_t PlaneWave_IBS::Repulsion(const ScalarFunction<double>& rho, double& Eh) const
{
    std::vector<rvec3_t> frac=UniformGrid(AutoGrid());
    UnitCell A=Recip().GetCell().MakeReciprocalCell();
    std::vector<double> field(frac.size());
    for (size_t q=0;q<frac.size();q++) field[q]=rho(A.ToCartesian(frac[q]));
    return Repulsion(ForwardDFTDiffSet(Gs(),frac,field), Eh);
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
    Eh *= 0.5*Volume();
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
    const std::vector<ivec3_t>& G=Gs();
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
            rg[G[i]-G[j]] += D(i,j);
    for (auto& kv : rg) kv.second /= Volume();
    return rg;
}

// Structure-factor assembly of a per-species radial form factor (the SAD seed density face): for each
// difference vector dm in the basis, rho(dm) = (1/Omega) Sum_atoms formFactor(Z,|B.dm|^2) e^{-i(B.dm).R}.
// Mirrors MakeLocalPotential, but it is a DENSITY: dm=0 is KEPT (= total charge / Omega), not dropped.
FourierMap PlaneWave_IBS::MakeFourierDensity(const Structure* atoms,
                          const std::function<double(int,double)>& formFactor) const
{
    const UnitCell& B=Recip().GetCell();
    FourierMap rho;
    size_t n=GetNumFunctions();
    const std::vector<ivec3_t>& G=Gs();
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            ivec3_t dm=G[i]-G[j];
            if (rho.find(dm)!=rho.end()) continue;     // one value per difference vector
            rvec3_t dG=B.ToCartesian(rvec3_t(dm));
            double  g2=dG*dG;
            dcmplx  acc(0.0);                           // (form factor) x (structure factor)
            for (Atom* a : *atoms) acc += formFactor(a->itsZ,g2)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
            rho[dm]=acc/Volume();
        }
    return rho;
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
    const std::vector<ivec3_t>& G=Gs();
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            ivec3_t dm=G[i]-G[j];
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
    return blazem::sum(f)*Volume()/double(f.size());
}

// Scalar integral integral f d3r over the cell: uniform-grid quadrature (weight Omega/Npts).
double PlaneWave_IBS::Integral(const ScalarFunction<double>& f) const
{
    std::vector<rvec3_t> frac=UniformGrid(AutoGrid());
    UnitCell A=Recip().GetCell().MakeReciprocalCell();
    double s=0.0;
    for (const rvec3_t& p : frac) s += f(A.ToCartesian(p));
    return s*Volume()/double(frac.size());
}

// Bare-Coulomb electron-nucleus attraction in reciprocal space: bare -Z/r is the bare-Coulomb local model.
//   <G|V|G'> = V(dG) = -(4 pi / Omega) Sum_a Z_a e^{-i dG.tau_a} / |dG|^2,   dG = G-G' != 0.
// The dG=0 term (the divergent G=0 Coulomb component) is dropped -- the conventional uniform background.
chmat_t PlaneWave_IBS::MakeNuclear(const Structure* cl) const
{
    return MakeLocalPotential(cl, Pseudopotential::BareCoulomb{});   // bare -Z/r is the bare-Coulomb local model
}

// V(dG) = (1/Omega) Sum_a v(Z_a,|dG|^2) e^{-i dG.tau_a}, dG=0 dropped (neutralising background).
// The result is Hermitian: V(-dG) = conj(V(dG)) since the structure factor conjugates under dG -> -dG
// (the real form factor v is even).
chmat_t PlaneWave_IBS::MakeLocalPotential(const Structure* cl, const Pseudopotential::LocalPotential& loc) const
{
    const UnitCell& B=Recip().GetCell();
    return MakePotential([&](const ivec3_t& dm)->dcmplx
    {
        if (dm.x==0 && dm.y==0 && dm.z==0) return dcmplx(0.0); // drop dG=0
        rvec3_t dG=B.ToCartesian(rvec3_t(dm));                  // dG = B.dm (Cartesian)
        double g2=dG*dG;
        dcmplx acc(0.0);                                        // (form factor) x (structure factor)
        for (Atom* a : *cl) acc += loc.FormFactor(a->itsZ,g2)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
        return acc/Volume();
    });
}

// V_NL(G,G') = (1/Omega) Sum_a e^{-i(G-G').tau_a} Sum_p (2l_p+1) P_{l_p}(cos gamma) betã_p(|k+G|) D_p
//              betã_p(|k+G'|),  gamma = angle(k+G, k+G').
// Per atom & projector & m this is rank-1: |beta> D <beta|.  Hermitian; real for atoms at the origin.
chmat_t PlaneWave_IBS::MakeSeparablePotential(const Structure* cl, const Pseudopotential::SeparablePotential& v) const
{
    const UnitCell& B=Recip().GetCell();
    size_t n=GetNumFunctions();
    const std::vector<ivec3_t>& G=Gs();
    Internal::KPlusG kg(B, kFrac(), G);             // Cartesian k+G, |k+G|, and cos(gamma)

    int maxL=0;                                     // highest projector channel present
    for (Atom* a : *cl)
        for (size_t p=0; p<v.NumProjectors(a->itsZ); p++)
            maxL=std::max(maxL, v.AngularMomentum(a->itsZ,p));

    chmat_t V=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            rvec3_t dG=B.ToCartesian(rvec3_t(G[i]-G[j]));
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
            V(i,j)=acc/Volume();
        }
    return V;
}

// The Band_FT_IBS factory seam: hand back a distinct auxiliary density-fit basis.  The density is
// cell-periodic, so ITS symmetry is Gamma (k=0), NOT this orbital block's k -- build a k=0 grid engine at
// the same Ecut and label it with a k=0 Bloch irrep.  (mp ignored -- a plane-wave fit basis is Ecut-based;
// a denser density grid is the future x2 tuning.)
BasisSet::cFIT_CD_ABS* PlaneWave_IBS::CreateCDFitBasisSet(const Structure*, const qcMesh::MeshParams&) const
{
    PW_Evaluator gamma(Recip(), rvec3_t(0.0,0.0,0.0), Ecut());
    return new PlaneWaveFit_IBS(gamma, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}

std::string PlaneWave_IBS::BasisSetID() const
{
    return Name()+PW_Evaluator::IDFragment();   // Name + "|k=..|Ecut=..|nG=.."
}

std::ostream& PlaneWave_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " plane waves, "
              << GetSymmetry();
}

} //namespace
