// File: src/BasisSet/Lattice/tests/PW_Fields.C  The plane-wave TEST FIXTURES that need no Hamiltonian:
// direct-grid oracles for the FFT/Poisson machinery (forward DFT, field overlaps, the Poisson solve), the
// hand-built rho~ map and a small cubic-cell fixture.  Shared by the PW basis unit tests (UTLattice_BS), the
// PW term tests (UTHamiltonian) and the PW SCF grid (IntegrationTests/PW) -- one copy, its own library
// (qcPW_TestFields), Hamiltonian-free so the Lattice test exe may link it.  (Moved out of the old
// PlaneWaveDFTUT.C on 2026-09-15, doc/TestSuitePlan.md phase 5.)
module;
#include <map>
#include <set>
#include <memory>
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cstdio>
export module qchem.Tests.PW_Fields;
export import qchem.BasisSet.PlaneWave.PlaneWave_IBS;
export import qchem.BasisSet.PlaneWave.Evaluators;   // PW_Grid_Evaluator -- a UNIT TEST may reach the internal
export import qchem.BasisSet.Lattice.BasisSet;   // Factory(Type::PW, lat, Ecut, loc, nl) -> Complex_BS*
export import qchem.ScalarFunction;                 // ScalarFunction<double> -- arg of the moved-here field oracles
export import qchem.Mesh;                           // qcMesh::MeshParams (Vee_Hartree's fit-basis factory arg; ignored)
export import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
export import qchem.Ewald;          // EwaldEnergy (ion-ion Madelung term -> physical total energy)
export import qchem.Types;          // dcmplx, ivec3_t, rvec_t, mat_t, chmat_t
export import qchem.Blaze;          // mat_t<dcmplx>
export import qchem.Math;           // Pi
export import qchem.Pseudopotential.GTH_Potentials;    // GetGTH (CP2K GTH/HGH database reader)
export import qchem.Fitting.FunctionFitter;                // Factory / ProjectedDensity_G / FunctionFitter_Density (item B)
export import qchem.BasisSet.G_FieldEvaluator;             // the grid-engine seam (GridPoints/RhoOnGrid/Integral) for the item-K probe
export import qchem.BasisSet.GMap;                             // ΔG_Map (the G-space coefficient map)
export import qchem.BasisSet.Orbital_DFT_IBS;                      // cFIT_CD_ABS (the ortho density-fit basis face)
export import qchem.Symmetry.Irrep;                        // Irrep
export import qchem.LASolver;                              // complex Hermitian eigensolver
export import qchem.Structure;                             // Molecule, Atom (the Si diamond basis)
export import qchem.Matrix3D;                              // Matrix3D<double> (the FCC cell matrix)
import qchem.BasisSet.Internal.BasisSetImp;         // BasisSetImp<dcmplx> (single-block BasisSet container)

export namespace qchem::tests::pw
{
using namespace qchem;
using BasisSet::PlaneWave::PlaneWave_IBS;
using BasisSet::PlaneWave::PW_Grid_Evaluator;   // internal grid evaluator (unit test may reach it directly)
using Pseudopotential::HGH_LocalPotential;
using Pseudopotential::HGH_SeparablePotential;
using Pseudopotential::GetGTH;
using Pseudopotential::GTH_PP;

// rho-tilde from a density matrix D via the basis's D-free Overlap3C tensor (the production path now that
// GetG_ERI3 is retired): Overlap3C keys on a Vxc fit basis (its grid is ignored -- the delta support is
// orbital-intrinsic).  Mirrors IrrepCD::GetFourierDensity(cFIT_SF_ABS).
ΔG_Map RhoTilde(const PlaneWave_IBS& pw, const chmat_t& D)
{
    std::unique_ptr<const qchem::BasisSet::cFIT_SF_ABS> vxcfb(pw.CreateVxcFitBasisSet(nullptr, qcMesh::MeshParams{}));
    return Contract(pw.Overlap3C(*vxcfb), D);
}
// Hartree matrix from a rho-tilde: <i|V_H|j> = 4pi/|G_i-G_j|^2 rho-tilde(G_i-G_j) (dm=0 dropped).  This is
// the retired Repulsion(ΔG_Map) route, inlined as a test cross-check (production assembles the same matrix
// from the density's Repulsion3C tensor via Vee_Hartree).
chmat_t HartreeFromRhoTilde(const PlaneWave_IBS& pw, const ΔG_Map& rt)
{
    return pw.OverlapMatrix([&](const ivec3_t& dm)->dcmplx
        { auto it=rt.find(dm); return it==rt.end()?dcmplx(0.0):pw.Recip().CoulombKernel(dm)*it->second; });
}

// A ScalarFunction<double> wrapping a lambda f(r) -- to hand real-space fields to the basis's
// high-level integral methods (Overlap/Repulsion/Integral).
struct FieldFn : public ScalarFunction<double>
{
    std::function<double(const rvec3_t&)> f;
    explicit FieldFn(std::function<double(const rvec3_t&)> g) : f(g) {}
    virtual double  operator()(const rvec3_t& r) const {return f(r);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}
};

// A v_xc(r) real field presented as a ProjectedScalar_R, so a scalar fitter can fit it (item K probe).
struct FieldFnR : public ScalarFunction<double>, public qchem::Fitting::ProjectedScalar_R
{
    std::function<double(const rvec3_t&)> f;
    explicit FieldFnR(std::function<double(const rvec3_t&)> g) : f(g) {}
    virtual double  operator()(const rvec3_t& r) const override {return f(r);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const override {return rvec3_t(0,0,0);}
    virtual const ScalarFunction<double>* GetScalarFunction() const override {return this;}
};

// --- Real-space DFT-integration oracles.  These lived on PlaneWave_IBS but were test-only cross-checks of
// the FFT/Poisson machinery (no library code called them), so they moved HERE as free functions over the
// basis's PUBLIC evaluator grid accessors (UniformGrid/AutoGrid/Gs/OverlapMatrix/Volume/Recip) -- no
// friendship needed, they never touched private data.

// Forward-transform a real field sampled on the fractional grid to its Fourier components over the difference
// set {m_i-m_j} (the only components the matrix <G_i|V|G_j> needs).
inline ΔG_Map ForwardDFTDiffSet(const std::vector<ivec3_t>& G, const std::vector<rvec3_t>& frac,
                                const std::vector<double>& field)
{
    ΔG_Map vt;
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

// The FFT/Poisson grid matching an orbital block: same (B,k,Ecut) -> same {G} -> same grid.  The orbital
// PlaneWave_IBS is grid-free now (the grid is a density/fit concern), so these direct-grid oracles build the
// density-grid evaluator themselves -- exactly the "test the evaluators directly" path.
inline qchem::BasisSet::PlaneWave::PW_Grid_Evaluator GridOf(const PlaneWave_IBS& pw)
{
    return qchem::BasisSet::PlaneWave::PW_Grid_Evaluator(pw.Recip(), pw.kFrac(), pw.Ecut());
}

// <i|f|j> weighted overlap for a real-space scalar field f: sample f on the (Cartesian) grid, forward-DFT to
// f-tilde(dm), assemble.  (The XC term's f(r)=v_xc(rho(r)) cross-check.)
inline chmat_t OverlapField(const PlaneWave_IBS& pw, const ScalarFunction<double>& f)
{
    std::vector<rvec3_t> frac=GridOf(pw).UniformGrid(pw.AutoGrid());
    UnitCell A=pw.Recip().GetCell().MakeReciprocalCell();        // direct cell (reciprocal of the reciprocal)
    std::vector<double> field(frac.size());
    for (size_t q=0;q<frac.size();q++) field[q]=f(A.ToCartesian(frac[q]));
    auto vt=ForwardDFTDiffSet(pw.Gs(),frac,field);
    return pw.OverlapMatrix([&vt](const ivec3_t& dm)->dcmplx
        { auto it=vt.find(dm); return it==vt.end()?dcmplx(0.0):it->second; });
}

// Coulomb repulsion matrix for a real-space density rho: sample on the grid, forward-DFT to rho-tilde, then
// the G-space Poisson solve V_H(dm)=4pi rho-tilde(dm)/|G|^2 assembled as <i|V_H|j>=V_H(G_i-G_j) (dm=0 dropped).
inline chmat_t RepulsionField(const PlaneWave_IBS& pw, const ScalarFunction<double>& rho)
{
    std::vector<rvec3_t> frac=GridOf(pw).UniformGrid(pw.AutoGrid());
    UnitCell A=pw.Recip().GetCell().MakeReciprocalCell();
    std::vector<double> field(frac.size());
    for (size_t q=0;q<frac.size();q++) field[q]=rho(A.ToCartesian(frac[q]));
    ΔG_Map rg=ForwardDFTDiffSet(pw.Gs(),frac,field);
    return pw.OverlapMatrix([&pw,&rg](const ivec3_t& dm)->dcmplx
    {
        auto it=rg.find(dm);
        return pw.Recip().CoulombKernel(dm)*(it==rg.end()?dcmplx(0.0):it->second);
    });
}

// integral f d3r over the cell: uniform-grid quadrature (weight Omega/Npts).
inline double IntegralField(const PlaneWave_IBS& pw, const ScalarFunction<double>& f)
{
    std::vector<rvec3_t> frac=GridOf(pw).UniformGrid(pw.AutoGrid());
    UnitCell A=pw.Recip().GetCell().MakeReciprocalCell();
    double s=0.0;
    for (const rvec3_t& p : frac) s += f(A.ToCartesian(p));
    return s*pw.Volume()/double(frac.size());
}

// Real-space grid values V -> matrix <i|V|j> = Vtilde(m_i-m_j).  Was PlaneWave_IBS::Overlap(rvec_t), now
// production-dead (the Vxc term assembles through the fit basis's seam); kept here as a test oracle: one
// ForwardFFT (the G_FieldEvaluator grid engine), then the orbital's OverlapMatrix lookup takes each
// reciprocal-index difference up in the grid.  (The orbital is both faces, so the PlaneWave_IBS supplies both.)
inline chmat_t OverlapOnGrid(const PlaneWave_IBS& pw, const rvec_t& V)
{
    PW_Grid_Evaluator grid=GridOf(pw);
    cvec_t Vt=grid.ForwardFFT(V);
    return pw.OverlapMatrix([&](const ivec3_t& dm)->dcmplx {return grid.GridCoeff(Vt, dm);});
}

// Order ivec3_t lexicographically so it can key the rho~ map.
struct IVecLess
{
    bool operator()(const ivec3_t& a, const ivec3_t& b) const
    {
        if (a.x!=b.x) return a.x<b.x;
        if (a.y!=b.y) return a.y<b.y;
        return a.z<b.z;
    }
};

//! rho~(dm): the periodic charge density's Fourier components, keyed by reciprocal-index difference dm.
typedef std::map<ivec3_t, dcmplx, IVecLess> RhoG;

//! rho~(dm) = (1/Omega) Sum_bands f_b Sum_{i,j} c_b(i) conj(c_b(j)) delta(m_i-m_j, dm).
//! U columns are the bands (c_b(i)=U(i,b)); f are the occupations; Omega is the direct-cell volume.
RhoG BuildDensity(const PlaneWave_IBS& pw, const mat_t<dcmplx>& U, const rvec_t& f, double Omega)
{
    size_t n=pw.GetNumFunctions();
    RhoG rho;
    for (size_t b=0; b<f.size(); b++)
    {
        if (f[b]==0.0) continue;
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                rho[pw.GetGIndex(i)-pw.GetGIndex(j)] += f[b]*U(i,b)*std::conj(U(j,b));
    }
    for (auto& kv : rho) kv.second/=Omega;
    return rho;
}

//! rho~(dm), or 0 if that component is absent (outside the difference set).
dcmplx RhoAt(const RhoG& rho, const ivec3_t& dm)
{
    auto it=rho.find(dm);
    return it==rho.end() ? dcmplx(0.0) : it->second;
}

// --- G-space Hartree:  V_H(r) solves nabla^2 V_H = -4 pi rho,  so V_H~(G) = 4 pi rho~(G)/|G|^2. -----
// The dG=0 component is dropped (neutralising background), as in MakeLocalPotential.  The matrix is
// then <G_i|V_H|G_j> = V_H~(m_i-m_j) via OverlapMatrix.

//! V_H~(dm) supplier for OverlapMatrix.  \a B is the RECIPROCAL cell (G = B dm).
std::function<dcmplx(const ivec3_t&)> HartreeVtilde(const RhoG& rho, const UnitCell& B)
{
    return [&rho,&B](const ivec3_t& dm)->dcmplx
    {
        if (dm==ivec3_t(0,0,0)) return dcmplx(0.0);
        rvec3_t G=B.ToCartesian(rvec3_t(dm));
        return 4*Pi*RhoAt(rho,dm)/(G*G);
    };
}

//! E_H = 1/2 integral rho V_H = (Omega/2) Sum_{G!=0} 4 pi |rho~(G)|^2 / |G|^2.
double HartreeEnergy(const RhoG& rho, const UnitCell& B, double Omega)
{
    double E=0.0;
    for (const auto& kv : rho)
    {
        if (kv.first==ivec3_t(0,0,0)) continue;
        rvec3_t G=B.ToCartesian(rvec3_t(kv.first));
        E += 4*Pi*std::norm(kv.second)/(G*G);
    }
    return 0.5*Omega*E;
}

// --- grid XC ----------------------------------------------------------------------------------
// LDA Vxc(rho(r)) is pointwise NONLINEAR, so it must be evaluated in real space: rho(G)->rho(r) on a
// uniform grid, apply the functional, then forward-transform Vxc(r)->Vxc~(G).  The grid is uniform
// (weight Omega/Npts): for the band-limited, cusp-free pseudo-density the trapezoidal rule is
// spectrally exact (see memory project_dft_upgrade_plan -- this is why PW codes use uniform FFT grids,
// not Becke quadrature).  Phases use dG.r = 2 pi dm.r_frac (since B = 2 pi A^{-T}), so the transform is
// a plain DFT independent of cell shape.  Direct DFT sums here; FFT is the later optimisation.

//! Uniform N1xN2xN3 grid of FRACTIONAL coordinates r = (i1/N1, i2/N2, i3/N3).
std::vector<rvec3_t> UniformGrid(const ivec3_t& Ng)
{
    std::vector<rvec3_t> g;
    g.reserve(size_t(Ng.x)*Ng.y*Ng.z);
    for (int i1=0;i1<Ng.x;i1++)
        for (int i2=0;i2<Ng.y;i2++)
            for (int i3=0;i3<Ng.z;i3++)
                g.push_back(rvec3_t(i1/double(Ng.x), i2/double(Ng.y), i3/double(Ng.z)));
    return g;
}

//! Inverse transform: rho(r) = Sum_dm rho~(dm) e^{i 2 pi dm.r_frac} (the imaginary part cancels: rho real).
double RhoOfR(const RhoG& spec, const rvec3_t& rf)
{
    dcmplx s(0.0);
    for (const auto& kv : spec)
    {
        double ph=2*Pi*(kv.first.x*rf.x + kv.first.y*rf.y + kv.first.z*rf.z);
        s += kv.second*dcmplx(std::cos(ph),std::sin(ph));
    }
    return std::real(s);
}

//! Forward transform: spectrum component f~(dm) = (1/Npts) Sum_grid f(r) e^{-i 2 pi dm.r_frac}.
dcmplx ForwardDFT(const std::vector<double>& f, const std::vector<rvec3_t>& rf, const ivec3_t& dm)
{
    dcmplx s(0.0);
    for (size_t q=0;q<f.size();q++)
    {
        double ph=-2*Pi*(dm.x*rf[q].x + dm.y*rf[q].y + dm.z*rf[q].z);
        s += f[q]*dcmplx(std::cos(ph),std::sin(ph));
    }
    return s/double(f.size());
}

// A small PW basis for the density tests: cubic cell, Gamma point.
struct PWFixture
{
    double           a=8.0, Ecut=4.0, Omega=a*a*a;
    UnitCell         cell{a};
    Lattice_3D       lat{cell, ivec3_t(1,1,1)};
    ReciprocalLattice recip{lat.Reciprocal()};   // owns the reciprocal cell B used by Hartree
    PlaneWave_IBS    pw{lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), Ecut};
    const UnitCell&  B() const {return recip.GetCell();}
};

//! Largest |component| over the G-set -- the difference set spans +/-2x this, so a grid with
//! N > 2*(2*maxComp) per axis resolves it without aliasing (=> band/direct energies agree exactly).
int MaxGComponent(const PlaneWave_IBS& pw)
{
    int m=0;
    for (size_t i=0;i<pw.GetNumFunctions();i++)
    {
        ivec3_t g=pw.GetGIndex(i);
        m=std::max(m, std::max({std::abs(g.x),std::abs(g.y),std::abs(g.z)}));
    }
    return m;
}

} // namespace qchem::tests::pw
