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
ΔG_Map
ForwardDFTDiffSet(const std::vector<ivec3_t>& G, const std::vector<rvec3_t>& frac,
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

// Coulomb repulsion matrix for a real-space density rho (TEST ORACLE): sample on the grid, forward-DFT to
// rho-tilde, then the G-space Poisson solve V_H(dm)=4pi rho-tilde(dm)/|G|^2 assembled as <i|V_H|j>=V_H(G_i-G_j)
// (dm=0 dropped).  The diagonal kernel is the reciprocal lattice's (Recip().CoulombKernel); MakePotential does
// the assembly.  (Production builds the Hartree matrix from the density's Repulsion3C tensor, not this route.)
chmat_t PlaneWave_IBS::Repulsion(const ScalarFunction<double>& rho) const
{
    std::vector<rvec3_t> frac=UniformGrid(AutoGrid());
    UnitCell A=Recip().GetCell().MakeReciprocalCell();
    std::vector<double> field(frac.size());
    for (size_t q=0;q<frac.size();q++) field[q]=rho(A.ToCartesian(frac[q]));
    ΔG_Map rg=ForwardDFTDiffSet(Gs(),frac,field);
    return MakePotential([this,&rg](const ivec3_t& dm)->dcmplx
    {
        auto it=rg.find(dm);
        return Recip().CoulombKernel(dm)*(it==rg.end()?dcmplx(0.0):it->second);
    });
}

// MakeRepulsion3C/MakeOverlap3C (the D-free {G} 3-centre tensor builds) moved to the shared grid engine
// (PW_Evaluator::Repulsion3CTensor/Overlap3CTensor); EPW_Orbital_DFT_IBS forwards the Band_FT_IBS virtuals to
// them.  MakeFourierDensity likewise moved (PW_Evaluator), so the SAD seed reaches it through its OWN fit basis.

// Scalar integral integral f d3r over the cell: uniform-grid quadrature (weight Omega/Npts).
double PlaneWave_IBS::Integral(const ScalarFunction<double>& f) const
{
    std::vector<rvec3_t> frac=UniformGrid(AutoGrid());
    UnitCell A=Recip().GetCell().MakeReciprocalCell();
    double s=0.0;
    for (const rvec3_t& p : frac) s += f(A.ToCartesian(p));
    return s*Volume()/double(frac.size());
}

// MakeNuclear (bare-Coulomb 1E block) moved to the evaluator (PW_Evaluator::NuclearMatrix), inherited via
// EPW_Orbital1E_IBS.

// The external LOCAL pseudopotential: forward to the shared structure-factor assembly on the grid engine
// (PW_Evaluator::LocalPotentialMatrix) with this model's form factor.  The Integrals_Pseudo capability stays
// here (it owns the LocalPotential model); the assembly loop lives ONCE on the evaluator (also drives MakeNuclear).
chmat_t PlaneWave_IBS::MakeLocalPotential(const Structure* cl, const Pseudopotential::LocalPotential& loc) const
{
    return LocalPotentialMatrix(cl, [&loc](int Z, double g2){ return loc.FormFactor(Z,g2); });
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

// The Band_FT_IBS factory seam: hand back a distinct auxiliary DENSITY-fit basis.  The density is
// cell-periodic, so ITS symmetry is Gamma (k=0), NOT this orbital block's k -- build a k=0 grid engine and
// label it with a k=0 Bloch irrep.  The electron density rho = psi*psi is EXACTLY band-limited to the
// difference set {G_i-G_j}: |G_i-G_j| <= 2 max|G|, i.e. 1/2|dG|^2 < 4 Ecut.  So the CD fit basis must cover
// 4x the orbital cutoff to represent rho -- honestly denser than the orbital {G} (never orbital==fit).  A
// GGA that wants a still-denser density grid raises it via mp.relCutoff (CP2K REL_CUTOFF); the 4x floor is
// the exact difference-set cover.  (The rho-tilde support is orbital-intrinsic today, so this is inert for the
// matrix -- the fitter holds it for the future denser-grid resampling -- but it makes the fit basis correct.)
BasisSet::cFIT_CD_ABS* PlaneWave_IBS::CreateCDFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    PW_Evaluator gamma(Recip(), rvec3_t(0.0,0.0,0.0), Ecut()*std::max(4.0, mp.relCutoff));
    return new PlaneWaveFit_IBS(gamma, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}

// The overlap-metric (Vxc) sibling: a DISTINCT fit-basis instance (same k=0 Gamma grid as the CD one when
// relCutoff matches; a GGA's denser Vxc grid diverges them via a larger GridCutoffFactor()).
// PlaneWaveFit_IBS IS-A cFIT_SF_ABS.
BasisSet::cFIT_SF_ABS* PlaneWave_IBS::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    PW_Evaluator gamma(Recip(), rvec3_t(0.0,0.0,0.0), Ecut()*mp.relCutoff);
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
