// File: BasisSet/PlaneWave/Imp/PlaneWave_IBS.C  Plane-wave irrep basis set implementation.
//
// Grid-geometry methods (op(r), overlap/kinetic, the {G} set, MakeOverlap, the FFT grid) live in the
// shared PW_Evaluator (this basis IS-A one, reached through the Lattice::*_IBS mixins).  What remains here is the
// orbital-only, atom/model-driven assembly: the density-driven G-space Hartree/XC route and the external
// pseudopotential.  Those read the shared grid data through the evaluator accessors (Gs(), Volume(),
// Recip(), kFrac(), GetGCartesian, MakeOverlap, the FFT-grid helpers).
module;
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

module qchem.BasisSet.PlaneWave.PlaneWave_IBS;
import qchem.BasisSet.PlaneWave.PlaneWaveFit_IBS;   // the auxiliary density-fit basis CreateCDFitBasisSet builds
import qchem.Symmetry.Factory;   // BlochFactory (the convenience ctor builds the Bloch irrep)
import qchem.Symmetry.Lattice_3D.BlochQN;   // Symmetry::Lattice_3D::Getk (prys k out of the abstract Bloch irrep)
import qchem.Structure;          // Atom (itsZ, itsR) + atom iteration for MakeNuclear
import qchem.Math;               // Pi, FourPi, cos, sin
import qchem.SpecialFunctions;   // LegendreP (the (2l+1)P_l angular factor)
import qchem.BasisSet.PlaneWave.Internal.KPlusG;     // KPlusG (Cartesian k+G, |k+G|, cos gamma)
import qchem.BasisSet.PlaneWave.Internal.GVectors;   // BuildGs (the enlarged sphere of MakeSetOverlap)
import qchem.Math.Angular;       // SphericalShell (the set's raw solid harmonics)
import qchem.Reporting;          // report::Timed (the +U set-overlap bucket)
import qchem.Blaze;
import qchem.Vector3D;           // dot product (operator*) + vector arithmetic

namespace qchem::BasisSet::PlaneWave
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

// The real-space DFT-integration oracles -- Overlap(f)/Repulsion(rho)/Integral(f) over a ScalarFunction --
// were test-only cross-checks; they moved to UnitTests/PlaneWaveDFTUT.C as free functions over the public
// evaluator grid accessors (UniformGrid/AutoGrid/Gs/MakeOverlap/Volume/Recip).  MakeRepulsion3C/MakeOverlap3C
// (the D-free {G} 3-centre builds) and MakeFourierDensity moved to the shared grid engine (PW_Evaluator).

// MakeNuclear (bare-Coulomb 1E block) moved to the evaluator (PW_Evaluator::NuclearMatrix), inherited via
// Lattice::Orbital_1E_IBS.

// A species LOCAL FIELD: forward to the shared structure-factor assembly on the grid engine
// (PW_Evaluator::LocalPotentialMatrix) with the field's form factor for the requested range.  The service
// face (Orbital_PP_IBS) stays here; the assembly loop lives ONCE on the evaluator (also drives MakeNuclear).
// A plane-wave basis assembles every range analytically in G-space (no grid), so the range is a trivial
// callback swap -- the CP2K split's efficiency payoff is GPW-side (it dodges the sharp-field local grid
// sweep for the long, deep-well part), but the uniform face lets PW_Hartree/PW_Pseudo drive either basis
// without a basis-type branch.
chmat_t PlaneWave_IBS::MakeSpeciesFieldMatrix(const Structure* cl, const SpeciesRadialField& f, FieldRange rng) const
{
    return LocalPotentialMatrix(cl, [&f,rng](int Z, double g2){ return f.ValueQ(Z,g2,rng); });
}

// V_NL(G,G') = (1/Omega) Sum_a e^{-i(G-G').tau_a} Sum_p (2l_p+1) P_{l_p}(cos gamma) betã_p(|k+G|) D_p
//              betã_p(|k+G'|),  gamma = angle(k+G, k+G').
// Per atom & projector & m this is rank-1: |beta> D <beta|.  Hermitian; real for atoms at the origin.
chmat_t PlaneWave_IBS::MakeProjectorMatrix(const Structure* cl, const SpeciesProjectorSet& v) const
{
    const UnitCell& B=Recip().GetCell();
    size_t n=GetNumFunctions();
    const std::vector<ivec3_t>& G=Gs();
    Internal::KPlusG kg(B, kFrac(), G);             // Cartesian k+G, |k+G|, and cos(gamma)

    int maxL=0;                                     // highest projector channel present
    for (Atom* a : *cl)
        for (size_t p=0; p<v.Count(a->itsZ); p++)
            maxL=std::max(maxL, v.L(a->itsZ,p));

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
                for (size_t p=0; p<v.Count(a->itsZ); p++)
                {
                    int l=v.L(a->itsZ,p);
                    s += (2*l+1)*P[l] * v.RadialQ(a->itsZ,p,kg.Norm(i))*v.Weight(a->itsZ,p)
                                       *v.RadialQ(a->itsZ,p,kg.Norm(j));
                }
                acc += s*std::exp(dcmplx(0.0,-(dG*a->itsR)));
            }
            V(i,j)=acc/Volume();
        }
    return V;
}

//------------------------------------------------------------------------------ the +U projector set (OpenWork PW-3)
namespace
{
// int over the unit sphere of x^a y^b z^c: 0 unless all even, else 2 G((a+1)/2) G((b+1)/2) G((c+1)/2) / G((a+b+c+3)/2).
double SphereMoment(int a, int b, int c)
{
    if (a%2 || b%2 || c%2) return 0.0;
    return 2.0*std::tgamma(0.5*(a+1))*std::tgamma(0.5*(b+1))*std::tgamma(0.5*(c+1))/std::tgamma(0.5*(a+b+c+3));
}
// Generalised Laguerre L_n^{(a)}(x) by the three-term recurrence.
double Laguerre(int n, double a, double x)
{
    double L0=1.0, L1=1.0+a-x;
    if (n==0) return L0;
    for (int k=1;k<n;k++) { const double L2=((2*k+1+a-x)*L1-(k+a)*L0)/(k+1); L0=L1; L1=L2; }
    return L1;
}
// One set MANIFOLD in the form its transform needs: the 2l+1 raw solid harmonics of the set's convention
// (SphericalShell(l)[l+m], each scaled POSITIVELY to unit norm on the sphere -- BasisSet.AtomCentredSet) sharing ONE
// radial and ONE phase, so a G vector costs one radial + one phase per manifold, not per function.
struct SetManifold
{
    rvec3_t centre; int l;
    std::vector<std::vector<qchem::Math::CartTerm>> harm;   // [m], scale folded in: harm[m](q) = q^l Y_m(q-hat)
    std::vector<Math::Gaussian> radial;
    // chi-tilde_{M,m}(q) for Cartesian q (no 1/sqrt(Omega)), m = -l..l into out[0..2l].
    void FT(const rvec3_t& q, dcmplx* out) const
    {
        const double q2=q*q;
        double rad=0.0;
        for (const auto& g : radial)
        {
            const double x=q2/(4.0*g.alpha);
            rad += g.c*std::sqrt(Pi)*std::tgamma(g.n+1.0)*Laguerre(g.n, l+0.5, x)*std::exp(-x)
                   /(std::pow(2.0,l+2)*std::pow(g.alpha, l+g.n+1.5));
        }
        static const dcmplx mi(0.0,-1.0);
        dcmplx il(1.0); for (int k=0;k<l;k++) il*=mi;   // (-i)^l
        const dcmplx common=FourPi*il*rad*std::exp(dcmplx(0.0,-(q*centre)));
        for (size_t m=0;m<harm.size();m++)
        {
            double poly=0.0;
            for (const auto& t : harm[m]) poly += t.c*std::pow(q.x,t.p.n)*std::pow(q.y,t.p.l)*std::pow(q.z,t.p.m);
            out[m]=common*poly;
        }
    }
};
std::vector<SetManifold> Manifolds(const BasisSet::AtomCentredSet& set)
{
    const auto* g=dynamic_cast<const BasisSet::AtomCentredSet_Gaussian*>(&set);
    if (!g) throw std::invalid_argument("PlaneWave_IBS: the +U set must carry the closed-Gaussian face (AtomCentredSet_Gaussian)");
    std::vector<SetManifold> f;
    for (size_t M=0;M<set.Count();M++)
    {
        const int l=set.L(M);
        const auto shell=qchem::Math::SphericalShell(l);
        SetManifold sm{set.Centre(M), l, {}, g->AsGaussians(M)};
        for (int m=-l;m<=l;m++)
        {
            auto harm=shell[size_t(l+m)];
            double n2=0.0;
            for (const auto& a : harm) for (const auto& b : harm) n2+=a.c*b.c*SphereMoment(a.p.n+b.p.n, a.p.l+b.p.l, a.p.m+b.p.m);
            for (auto& t : harm) t.c/=std::sqrt(n2);
            sm.harm.push_back(std::move(harm));
        }
        f.push_back(std::move(sm));
    }
    return f;
}
// All N set functions at one q, manifold by manifold, m ascending (the set's column order).
void TransformAll(const std::vector<SetManifold>& f, const rvec3_t& q, dcmplx* out)
{
    for (const auto& M : f) { M.FT(q, out); out+=M.harm.size(); }
}
size_t NumFunctions(const std::vector<SetManifold>& f) { size_t n=0; for (const auto& M : f) n+=M.harm.size(); return n; }
// F(G,a) = chi-tilde_a(k+G) over a G list.
mat_t<dcmplx> Transforms(const std::vector<SetManifold>& f, const UnitCell& B, const rvec3_t& k, const std::vector<ivec3_t>& G)
{
    const size_t N=NumFunctions(f);
    mat_t<dcmplx> F(G.size(), N);
    std::vector<dcmplx> v(N);
    for (size_t i=0;i<G.size();i++)
    {
        TransformAll(f, B.ToCartesian(k+rvec3_t(G[i])), v.data());
        for (size_t a=0;a<N;a++) F(i,a)=v[a];
    }
    return F;
}
} //anon

mat_t<dcmplx> PlaneWave_IBS::MakeAtomCentredOverlap(const AtomCentredSet& set) const
{
    mat_t<dcmplx> B=Transforms(Manifolds(set), Recip().GetCell(), kFrac(), Gs());
    B*=1.0/std::sqrt(Volume());
    return B;
}

chmat_t PlaneWave_IBS::MakeSetOverlap(const AtomCentredSet& set) const
{
    // ACCUMULATED, never stored: the pseudo-atomic radials carry sharp pool exponents, so 1e-12 needs ~1e7 G vectors,
    // and an N_G x N transform matrix (+ its adjoint for F^dag F) was 7.4 GB at 12 Ha (2026-10-25, scripts/bigalloc).
    // The sphere grows in SHELLS (cutoff doubled); each pass adds only the G new to it, and stops when a whole shell
    // moves O by < 1e-12.
    qchem::report::Timed timed("setup: +U set overlap (Parseval, accumulated by shell)");
    const auto f=Manifolds(set);
    const size_t N=NumFunctions(f);
    const UnitCell& B=Recip().GetCell();
    mat_t<dcmplx> O(N,N,dcmplx(0.0));
    std::vector<dcmplx> v(N);
    auto addShell=[&](double Elo, double Ehi)   // the G with Elo <= 1/2|k+G|^2 < Ehi; returns max |dO|
    {
        const std::vector<ivec3_t> G=Internal::BuildGs(Recip(), kFrac(), Ehi);
        if (G.size()>40000000)
            throw std::runtime_error("PlaneWave_IBS::MakeSetOverlap: the +U set's overlap needs > 4e7 G vectors (Ecut "
                                     +std::to_string(Ehi)+" Ha) and has not converged -- a radial too sharp for a Parseval sum");
        mat_t<dcmplx> dO(N,N,dcmplx(0.0));
        for (const ivec3_t& g : G)
        {
            const rvec3_t q=B.ToCartesian(kFrac()+rvec3_t(g));
            if (0.5*(q*q)<Elo) continue;
            TransformAll(f, q, v.data());
            for (size_t a=0;a<N;a++) { const dcmplx ca=std::conj(v[a]); for (size_t b=a;b<N;b++) dO(a,b)+=ca*v[b]; }
        }
        double d=0.0;
        for (size_t a=0;a<N;a++) for (size_t b=a;b<N;b++) { O(a,b)+=dO(a,b)/Volume(); d=std::max(d, std::abs(dO(a,b))/Volume()); }
        return d;
    };
    double E=std::max(Ecut(), 10.0);
    addShell(-1.0, E);
    while (addShell(E, 2.0*E)>=1e-12) E*=2.0;
    chmat_t H=blazem::zeroH<dcmplx>(N);
    for (size_t a=0;a<N;a++) for (size_t b=a;b<N;b++) H(a,b)= a==b ? dcmplx(std::real(O(a,b)),0.0) : O(a,b);
    return H;
}

// The Orbital_DFT_IBS<dcmplx> factory seam: hand back a distinct auxiliary DENSITY-fit basis.  The density is
// cell-periodic, so ITS symmetry is Gamma (k=0), NOT this orbital block's k -- build a k=0 grid engine and
// label it with a k=0 Bloch irrep.  The electron density rho = psi*psi is EXACTLY band-limited to the
// difference set {G_i-G_j}: |G_i-G_j| <= 2 max|G|, i.e. 1/2|dG|^2 < 4 Ecut.  So the CD fit basis must cover
// 4x the orbital cutoff to represent rho -- honestly denser than the orbital {G} (never orbital==fit).  A
// GGA that wants a still-denser density grid raises it via mp.relCutoff (CP2K REL_CUTOFF); the 4x floor is
// the exact difference-set cover.  (The rho-tilde support is orbital-intrinsic today, so this is inert for the
// matrix -- the fitter holds it for the future denser-grid resampling -- but it makes the fit basis correct.)
BasisSet::cFIT_CD_ABS* PlaneWave_IBS::CreateCDFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    PW_Grid_Evaluator gamma(Recip(), rvec3_t(0.0,0.0,0.0), Ecut()*std::max(4.0, mp.relCutoff));
    return new PlaneWaveFit_IBS(gamma, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)), "densityFit");
}

// The overlap-metric (Vxc) sibling: a DISTINCT fit-basis instance (same k=0 Gamma grid as the CD one when
// relCutoff matches; a GGA's denser Vxc grid diverges them via a larger GridCutoffFactor()).
// PlaneWaveFit_IBS IS-A cFIT_SF_ABS.
BasisSet::cFIT_SF_ABS* PlaneWave_IBS::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    // A BLOCK builds only its OWN fitted representation (see the GPW sibling): delta-vs-fitted is the
    // whole-set tBasisSet factory's decision.
    PW_Grid_Evaluator gamma(Recip(), rvec3_t(0.0,0.0,0.0), Ecut()*mp.relCutoff);
    return new PlaneWaveFit_IBS(gamma, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)), "vxcFitGrid");
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
