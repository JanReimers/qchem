// File: BasisSet/Lattice_3D/Evaluators/PW/Imp/Evaluator.C  PW_Evaluator implementation.
module;
#include <cassert>
#include <algorithm>
#include <cstdlib>    // getenv (GPW_FIELD_SPECTRUM)
#include <iostream>
#include <complex>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.Evaluators.PW;
import qchem.Math;               // Pi, sqrt, cos, sin, Cube
import qchem.BasisSet.Lattice_3D.Internal.GVectors;   // BuildGs
import qchem.FFT;                // NextPow2 (the XC grid geometry)
import qchem.Blaze;              // blazem::zeroH
import qchem.Vector3D;           // dot product (operator*) + vector arithmetic

namespace qchem::BasisSet::Lattice_3D
{

PW_Evaluator::PW_Evaluator(const ReciprocalLattice& recip, const rvec3_t& k, double Ecut,
                           RasterPolicy raster)
    : itsRecip(recip)
    , itsk(k)
    , itsEcut(Ecut)
    // |det B| = (2 pi)^3 / |det A|, so the direct cell volume V = (2 pi)^3 / V_recip.
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
    , itsRaster(raster)
{
    itsG = Internal::BuildGs(itsRecip, itsk, Ecut);   // { G : 1/2|k+G|^2 < Ecut }
}

rvec3_t PW_Evaluator::GetGCartesian(const ivec3_t& m) const
{
    return itsRecip.GetCell().ToCartesian(rvec3_t(m)); // B m
}

// Grid divisions resolving the difference set without aliasing: N > 2*(2*maxComp).  Cached (depends on itsG).
ivec3_t PW_Evaluator::AutoGrid() const
{
    if (itsAutoGrid.x==0)   // not yet computed
    {
        int m=0;
        for (const ivec3_t& g : itsG)
        {
            int ax=g.x<0?-g.x:g.x, ay=g.y<0?-g.y:g.y, az=g.z<0?-g.z:g.z;
            m=std::max(m, std::max(ax, std::max(ay,az)));
        }
        // AliasFree: 4m+1 resolves the difference set {G-G'} exactly; BallOnly: 2m+1 resolves the ball
        // alone, product fold-back accepted as Ecut-controlled discretization error (doc/GPWPlan 0.5(a)).
        int nn = itsRaster==RasterPolicy::AliasFree ? 4*m+1 : 2*m+1;
        itsAutoGrid=ivec3_t(nn,nn,nn);
    }
    return itsAutoGrid;
}

// AutoGrid divisions rounded up to 5-SMOOTH numbers (2^a 3^b 5^c) -- the mixed-radix FFT menu
// (qchem.FFT dispatches non-pow2 N to PocketFFT; doc/GPWPlan.md mixed-radix increment 2026-07-22).
// A larger grid still resolves the difference set without aliasing, so padding is at worst slightly
// more accurate; the FINER menu is the raster-volume lever (pow2 padding cost NaF 128^3 where the
// 5-smooth menu admits 75^3-class rasters -- every N^3-scaled cost shrinks: boxes, streams, sweeps).
// Flipping this policy moved every grid-derived anchor (the one-time re-pin, same increment).  Cached.
ivec3_t PW_Evaluator::FFTGrid() const
{
    if (itsFFTGrid.x==0)   // not yet computed
    {
        ivec3_t a=AutoGrid();
        itsFFTGrid=ivec3_t(int(qchem::FFT::Next5Smooth(a.x)), int(qchem::FFT::Next5Smooth(a.y)), int(qchem::FFT::Next5Smooth(a.z)));
    }
    return itsFFTGrid;
}

// The fitted field's coefficients over this basis's own {G} (a GridCoeff gather) -- for op(r) evaluation.
// A DENSITY/FIT operation (iterates {G} + the FFT grid), so it lives on PW_Grid_Evaluator, reaching {G} via
// the base accessor Gs() and the per-G lookup via the held grid engine.
// GPW_FIELD_SPECTRUM=1: HOW BIG DOES THIS {G} BALL ACTUALLY NEED TO BE?
//
// The Vxc fit basis is built at {G}_vxc = relCutoff * {G}_rho and relCutoff is 1.0 for every functional
// today (GridCutoffFactor defaults to 1), so v_xc is expanded over the FULL density ball -- 24k vectors on
// MnO.  That matters beyond memory: a plane-wave fit on a NON-RASTER mesh is no longer discretely
// orthonormal, so it needs the overlap metric S = <c_a|c_b>, and |{G}|^2 at 24k is not a matrix anyone
// inverts.  At 1/3 the RADIUS it is 1/27 the volume -- ~900 vectors, S ~ 900^2, one Cholesky per geometry.
//
// So the question is where the coefficients have actually decayed, and it is answerable for free from the
// fit that just ran.  Reported as the CUMULATIVE FRACTION OF sum|c|^2 inside a given radius fraction: read
// off the radius that captures 1-eps and that is the ball the fit needs.  Physically v_xc ~ rho^(1/3) for
// Dirac exchange, and the cube root flattens rho's cusps, so its spectrum should decay much faster than the
// density's -- this measures whether it does rather than assuming it.
void ReportFieldSpectrum(const ΔG_Map& m, const ReciprocalLattice& recip)
{
    static const bool on=std::getenv("GPW_FIELD_SPECTRUM")!=nullptr;
    if (!on || m.empty()) return;
    double gmax=0.0, total=0.0;
    for (const auto& [dm,c] : m) { gmax=std::max(gmax, recip.GetGLength(dm)); total+=std::norm(dcmplx(c)); }
    if (gmax<=0.0 || total<=0.0) return;
    constexpr size_t kNBin=10;
    double binPow[kNBin]={0,0,0,0,0,0,0,0,0,0}; size_t binN[kNBin]={0,0,0,0,0,0,0,0,0,0};
    for (const auto& [dm,c] : m)
    {
        size_t b=size_t((recip.GetGLength(dm)/gmax)*double(kNBin)); if (b>=kNBin) b=kNBin-1;
        binPow[b]+=std::norm(dcmplx(c)); binN[b]++;
    }
    std::cout<<"[field spectrum] |{G}|="<<m.size()<<" Gmax="<<gmax
             <<"  cumulative fraction of sum|c|^2 inside radius x*Gmax:";
    double cum=0.0; size_t cumN=0;
    for (size_t b=0;b<kNBin;b++)
    {
        cum+=binPow[b]/total; cumN+=binN[b];
        std::cout<<"  "<<(b+1)*10<<"%:"<<cum<<"("<<cumN<<"G)";
    }
    std::cout<<std::endl;
}

ΔG_Map PW_Grid_Evaluator::FieldCoeffs(const cvec_t& Vt) const
{
    ΔG_Map m;
    for (const ivec3_t& G : Gs()) m[G]=itsGrid->GridCoeff(Vt, G);
    ReportFieldSpectrum(m, Recip());
    return m;
}

// The ADJOINT of EvalField (see the face): c(G) = (1/Omega) Sum_g w_g v(r_g) e^{-iG.r_g}, over THIS basis's
// own {G} -- one weight-vector contraction per G, so it takes any point set and any weights.  O(|{G}| x npts)
// by construction: the FFT that does all of these at once needs the points to BE a raster, which is exactly
// the assumption this entry point exists to drop.  TRUNCATING, not aliasing (see the \warning on the face).
ΔG_Map PW_Grid_Evaluator::ProjectField(const rvec3vec_t& pts, const rvec_t& w, const rvec_t& vals) const
{
    assert(pts.size()==vals.size() && pts.size()==w.size() &&
           "PW_Grid_Evaluator::ProjectField: points, weights and values must agree in length");
    const UnitCell& B=Recip().GetCell();
    const double invOmega=1.0/Volume();
    ΔG_Map m;
    for (const ivec3_t& dm : Gs())
    {
        const rvec3_t G=B.ToCartesian(rvec3_t(double(dm.x),double(dm.y),double(dm.z)));
        double re=0.0, im=0.0;
        for (size_t g=0; g<pts.size(); ++g)
        {
            const double ph=G*pts[g], wv=w[g]*vals[g];
            re+=wv*cos(ph); im-=wv*sin(ph);          // e^{-iG.r}
        }
        m[dm]=invOmega*dcmplx(re,im);
    }
    return m;
}

// Analytic structure-factor density over this engine's own {G}: rho-tilde(G) = (1/Omega) Sum_atoms
// formFactor(Z,|B.G|^2) e^{-i(B.G).R}.  formFactor is the atomic density's 1-D radial Fourier transform, so
// this is ANALYTIC per G -- no 3-D grid, no aliasing of the peaked density (unlike sample+FFT).  A DENSITY:
// G=0 is KEPT (= total charge / Omega).  The SAD seed calls this on its OWN fit basis (never the orbital one).
// Grid-free ({G}+B+Omega from the base) but a density-side operation, so it lives on PW_Grid_Evaluator.
ΔG_Map PW_Grid_Evaluator::MakeFourierDensity(const Structure* atoms,
                          const std::function<double(int,double)>& formFactor) const
{
    const UnitCell& B=Recip().GetCell();
    ΔG_Map rho;
    for (const ivec3_t& dm : Gs())                     // this basis's own {G} (the density's Fourier support)
    {
        rvec3_t dG=B.ToCartesian(rvec3_t(dm));
        double  g2=dG*dG;
        dcmplx  acc(0.0);                               // (form factor) x (structure factor)
        for (Atom* a : *atoms) acc += formFactor(a->itsZ,g2)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
        rho[dm]=acc/Volume();
    }
    return rho;
}

// Reusable local-potential assembly: <G|V|G'> = (1/Omega) Sum_a f(Z_a,|dG|^2) e^{-i dG.tau_a}, dG=0 dropped
// (the neutralising background).  Hermitian: V(-dG)=conj(V(dG)) since the structure factor conjugates and the
// (real) form factor is even.  MakeNuclear is this with the bare-Coulomb f; a pseudopotential local term reuses
// it with its own form factor -- so the structure-factor loop lives ONCE, here on the grid engine.
chmat_t PW_Evaluator::LocalPotentialMatrix(const Structure* cl,
                          const std::function<double(int,double)>& formFactor) const
{
    const UnitCell& B=itsRecip.GetCell();
    return OverlapMatrix([&](const ivec3_t& dm)->dcmplx
    {
        if (dm.x==0 && dm.y==0 && dm.z==0) return dcmplx(0.0);   // drop dG=0
        rvec3_t dG=B.ToCartesian(rvec3_t(dm));                    // dG = B.dm (Cartesian)
        double  g2=dG*dG;
        dcmplx  acc(0.0);                                         // (form factor) x (structure factor)
        for (Atom* a : *cl) acc += formFactor(a->itsZ,g2)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
        return acc/Volume();
    });
}

// Bare-Coulomb electron-nucleus attraction: the local potential with f(Z,g2) = -4 pi Z / g2, i.e.
// <G|V|G'> = -(4 pi/Omega) Sum_a Z_a e^{-i dG.tau_a}/|dG|^2, dG=0 dropped.  The 1E nuclear block.
chmat_t PW_Evaluator::NuclearMatrix(const Structure* cl) const
{
    return LocalPotentialMatrix(cl, [](int Z, double g2){ return -FourPi*Z/g2; });
}

namespace
{
// The delta support: one column per difference dm=G_i-G_j, listing the (i,j) pairs that hit it, row-major
// (i outer, j inner) so a per-column left-fold reproduces the density contraction exactly.
void BuildProjector3Columns(const std::vector<ivec3_t>& G, std::vector<Projector3<dcmplx>::Column>& cols)
{
    size_t n=G.size();
    cols.clear();
    std::map<ivec3_t,int,IVec3Less> colOf;          // dm -> column index (build-time lookup only)
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            ivec3_t dm=G[i]-G[j];
            auto it=colOf.find(dm);
            int c;
            if (it==colOf.end()) { c=int(cols.size()); colOf[dm]=c; cols.push_back({dm,{}}); }
            else                   c=it->second;
            cols[c].pairs.push_back({int(i),int(j)});
        }
}
} //anon

// Coulomb 3-centre tensor over this engine's {G}: delta support + diagonal Poisson kernel 4pi/|G_c|^2
// (dm=0 -> 0, the dropped G=0 background).  A density contracts D against it (Contract) for V_H.
Projector3<dcmplx> PW_Evaluator::Repulsion3CTensor() const
{
    Projector3<dcmplx> g;
    BuildProjector3Columns(itsG, g.columns);
    g.kernel.resize(g.columns.size());
    for (size_t c=0;c<g.columns.size();c++)
        g.kernel[c]=itsRecip.CoulombKernel(g.columns[c].dm);   // diagonal Poisson kernel (dm=0 -> 0)
    g.volume=Volume();
    g.applyAdjoint=AdjointLookup();   // <i|f|j> = f(G_i-G_j) (overlap metric; the kernel is forward-only)
    return g;
}

// Overlap 3-centre tensor: delta support, EMPTY kernel (overlap metric).
Projector3<dcmplx> PW_Evaluator::Overlap3CTensor() const
{
    Projector3<dcmplx> g;
    BuildProjector3Columns(itsG, g.columns);
    g.kernel.clear();                // EMPTY => overlap metric
    g.volume=Volume();
    g.applyAdjoint=AdjointLookup();
    return g;
}

// The Projector3<dcmplx> BACKWARD realization for plane waves: <i|f|j> = f(m_i-m_j) (== OverlapMatrix), the Fourier lookup.
// Self-contained (captures the orbital {G} by value), so the closure outlives the evaluator in the tensor cache.
std::function<chmat_t(const std::function<dcmplx(const ivec3_t&)>&)> PW_Evaluator::AdjointLookup() const
{
    return [G=itsG](const std::function<dcmplx(const ivec3_t&)>& f) -> chmat_t
    {
        const size_t n=G.size();
        chmat_t V=blazem::zeroH<dcmplx>(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) V(i,j)=f(G[i]-G[j]);
        return V;
    };
}

// Plane waves are orthonormal over the cell: <G|G'> = delta_{GG'}.
chmat_t PW_Evaluator::OverlapMatrix() const
{
    size_t n=size();
    chmat_t S=blazem::zeroH<dcmplx>(n);   // hmat_t(n) does NOT zero the off-diagonals
    for (size_t i=0; i<n; i++) S(i,i)=1.0;
    return S;
}

// <p^2> = <-nabla^2> building block (NO 1/2 -- the Hamiltonian applies it).  For a plane wave
// -nabla^2 e^{i(k+G).r} = |k+G|^2 e^{i(k+G).r}, so the matrix is diagonal in |k+G|^2.
chmat_t PW_Evaluator::KineticMatrix() const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=size();
    chmat_t S=blazem::zeroH<dcmplx>(n);   // off-diagonals are exactly zero
    for (size_t i=0; i<n; i++)
    {
        double kG=B.GetDistance(itsk+itsG[i]); // |k+G|
        S(i,i)=kG*kG;
    }
    return S;
}

// <G|V|G'> = Vtilde(m(G) - m(G')).  Fill the upper triangle; HermitianMatrix mirrors the conjugate.
chmat_t PW_Evaluator::OverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    size_t n=size();
    chmat_t V=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
            V(i,j)=Vtilde(itsG[i]-itsG[j]);
    return V;
}

cvec_t PW_Evaluator::Eval(const rvec3_t& r) const
{
    size_t n=size();
    double invSqrtV=1.0/sqrt(itsVolume);
    cvec_t v(n);
    for (size_t i=0; i<n; i++)
    {
        double phase=GetGCartesian(itsG[i])*r + itsRecip.GetCell().ToCartesian(itsk)*r; // (k+G).r
        v[i]=dcmplx(cos(phase),sin(phase))*invSqrtV;
    }
    return v;
}

cvec3vec_t PW_Evaluator::EvalGradient(const rvec3_t& r) const
{
    // grad e^{i(k+G).r}/sqrt(V) = i(k+G) e^{i(k+G).r}/sqrt(V).
    const dcmplx im(0.0,1.0);
    size_t n=size();
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

std::string PW_Evaluator::IDFragment() const
{
    // The cache key must satisfy the DBCacheClient contract: equal physics <=> equal string.  k/Ecut/nG
    // alone are NOT enough -- every G-space matrix (kinetic 1/2|k+G|^2, the Coulomb kernel 4pi/|G|^2, the
    // 1/Omega normalisation) depends on the reciprocal CELL METRIC B, and two differently-shaped cells can
    // share an nG.  So append B's three Cartesian columns (B.e_i = ToCartesian of the unit fractional axes),
    // which pin |G|^2 and Omega exactly.  (Without this the process-wide cache could hand a same-nG neighbour
    // the wrong kernel -- the per-instance member cache this replaces could never collide.)
    const UnitCell& B=itsRecip.GetCell();
    auto v=[](const rvec3_t& c){ return std::to_string(c.x)+","+std::to_string(c.y)+","+std::to_string(c.z); };
    return std::string("|k=")+std::to_string(itsk.x)+","+std::to_string(itsk.y)+","+std::to_string(itsk.z)
                 +"|Ecut="+std::to_string(itsEcut)
                 +"|nG="+std::to_string(itsG.size())
                 +"|B="+v(B.ToCartesian(rvec3_t(1,0,0)))+";"+v(B.ToCartesian(rvec3_t(0,1,0)))+";"+v(B.ToCartesian(rvec3_t(0,0,1)));
}

} //namespace
