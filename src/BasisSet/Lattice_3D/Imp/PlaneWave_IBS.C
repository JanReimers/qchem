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

// The real-space DFT-integration oracles -- Overlap(f)/Repulsion(rho)/Integral(f) over a ScalarFunction --
// were test-only cross-checks; they moved to UnitTests/PlaneWaveDFTUT.C as free functions over the public
// evaluator grid accessors (UniformGrid/AutoGrid/Gs/MakePotential/Volume/Recip).  MakeRepulsion3C/MakeOverlap3C
// (the D-free {G} 3-centre builds) and MakeFourierDensity moved to the shared grid engine (PW_Evaluator).

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
    PW_Grid_Evaluator gamma(Recip(), rvec3_t(0.0,0.0,0.0), Ecut()*std::max(4.0, mp.relCutoff));
    return new PlaneWaveFit_IBS(gamma, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}

// The overlap-metric (Vxc) sibling: a DISTINCT fit-basis instance (same k=0 Gamma grid as the CD one when
// relCutoff matches; a GGA's denser Vxc grid diverges them via a larger GridCutoffFactor()).
// PlaneWaveFit_IBS IS-A cFIT_SF_ABS.
BasisSet::cFIT_SF_ABS* PlaneWave_IBS::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    PW_Grid_Evaluator gamma(Recip(), rvec3_t(0.0,0.0,0.0), Ecut()*mp.relCutoff);
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
