// File: IntegrationTests/GPW/MnO.C  MnO rocksalt AFM-II -- the magnetic transition-metal cell (rhombohedral 2-f.u., Mn +/-m, O); SEED-LEVEL gates only, the campaign run is CLIapps/gpwprobe mno.
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(<Basis>_<Material>, <k>_[Grid]_[Fit]_[Sym]_[Spin]_[Occ]_[Machinery]_[Ansatz]_[Seed]_<Claim>)
// -- axis tokens in fixed order, the facade's defaults elided (k always named); the claim is CP2K (oracle anchor),
// Anchor (did-E-move), eq<Token> (a ONE-axis twin: this point equals the same point with that axis moved) or a
// property verb.  Tests are laid out in axis order.  `scripts/testgrid` renders the coverage table from these names.
//
//   GPW_MnO.Γ_Pol_SeedMirror
//   GPW_MnO.Γ_Becke_Pol_SeedVxcMirror
//   GPW_MnO.Γ_Shub_Pol_SeedDecoration

#include "gtest/gtest.h"
#include <memory>
#include <vector>
#include <cmath>
#include <cstdlib>   // std::getenv/std::atof (the NaF mixing-tuning env knobs)
#include <complex>
#include <cstdio>
#include <fstream>   // /proc/self/statm (the RSS breadcrumb bisect)
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <string>
#include <iomanip>      // setprecision (the order-parameter trajectory line)

import qchem.Structure;                          // Molecule, Atom
import qchem.UnitCell;                           // UnitCell, FCCUnitCell
import qchem.Matrix3D;                           // Matrix3D<double> (the rhombohedral AFM-II cell matrix)
import qchem.Lattice_3D;                         // Lattice_3D
import qchem.BasisSet;                           // Complex_BS, Real_BS
import qchem.BasisSet.Orbital_1E_IBS;            // Complex_OIBS (the overlap-spectrum diagnostic)
import qchem.Blaze;                              // blazem::eigen, blaze::min/max (overlap spectrum)
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;   // PlaneWave_IBS (the seed's CD fit basis)
import qchem.BasisSet.Lattice.BasisSet;       // GPWFactory (the GPW basis container)
import qchem.BasisSet.Gaussian.Point.Factory;          // Gaussian::Factory, BasisSetData/Engine/Angular
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;  // MakeSphericalLatticeView (GPW_SPHERICAL=1)
import qchem.Hamiltonian.Factory;                 // the PUBLIC solid front door (Step 4): cHamiltonian* Factory(...)
import qchem.Outcome;                           // Outcome<Converged,SCFFailure> -- the facade's result
import qchem.RunPolicy;                         // ReresolveRunPolicy() -- the declared-deviation A/B hatch (N5)
import qchem.SolidCalculation;                    // the NAMED periodic facade (Step 4 3/3)
import qchem.Tests.GPW_Harness;                   // THE HARNESS (IntegrationTests/GPW/Harness.C): Materials cells, gates, recipes, the XC probes
import qchem.Materials;                           // Materials::Get -- the cells come from src/Calculation/Data/materials.json (row MD)
import qchem.Hamiltonian.Internal.Hamiltonians;  // Ham_PW_DFT direct ctors (the bespoke probes below still use them)
import qchem.Hamiltonian.Internal.PWTerms;        // ReportGridCharge(); Vxc_Quadrature + the two DensitySampler strategies
import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
                                                  // since 2026-09-08; .Internal. modules are
                                                  // never re-exported, so name it directly)
import qchem.BasisSet.DeltaFit_IBS;              // DeltaFit_IBS -- the delta basis the singles strategy runs on
import qchem.BasisSet.G_FieldEvaluator;           // G_RasterTransform -- the uniform probe's own point count
import qchem.Mesh.Angular;                        // MakeAngular (the rotated-Lebedev bond-angle probe)
import qchem.Hamiltonian.Internal.ExFunctional;   // ExFunctional (the LDA functional face the XC terms hold)
import qchem.Hamiltonian.Internal.SlaterExchange; // SlaterExchange (Dirac exchange, for the Becke XC gate)
import qchem.Hamiltonian.Internal.VWN_Correlation;// VWN_Correlation (VWN5, for the Becke XC gate)
import qchem.Mesh;                                // qcMesh::MeshParams / UnitCellKind (the Becke XC quadrature)
import qchem.Mesh.XCPolicy;                       // BeckeXCParams / ResolveXCMesh / XCMeshSharpness (the grid policy)
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;      // Gaussian::LatticeSum1E::MaxExponent (alpha_max, for the selector)
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH -> HGH local PP (alpha_pp = 1/2r_loc^2, for the selector)
import qchem.PeriodicTable;                       // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.SCFIterator;                        // cSCFIterator, SCFParams
import qchem.SCFParams;                          // SCFParams
import qchem.ElectronConfiguration.Crystal;      // Crystal_EC (single-k Bloch occupation)
import qchem.ChargeDensity.Seed;                 // SeedStrategy
import qchem.SCFAccelerator.Factory;              // the PUBLIC complex accelerator door (Step 4)
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS; // SCFAcceleratorDIIS (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;  // SCFAcceleratorGDM (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorLadder; // SCFAcceleratorLadder (DIIS -> GDM chain)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // SCFAcceleratorNull (NaF: pure damped Kerker)
import qchem.WaveFunction;                       // cWaveFunction (the converged state)
import qchem.Energy;                             // EnergyBreakdown
import qchem.Symmetry.Irrep;                     // Irrep
import qchem.Reporting;                          // report:: -- bracket the GPW run so grids/basis sections land
import qchem.Symmetry.Spin;                      // Spin
import qchem.Symmetry.Factory;                   // BlochFactory (build a k-block with a fractional MP shift)
import qchem.LASolver;                           // qchem::Ortho (Cholesky | Eigen | SVD -- basis orthogonalisation)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;         // GPW_IBS (build a concrete block for the collocation diagnostic)
import qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator;  // GPW_Evaluator (Overlap3CTensor -- the collocation tensor)
import qchem.BasisSet.GMap;              // Projector3<dcmplx> (the collocation weight tensor); SymmetryDefects (§3 diagnostic)
import qchem.ChargeDensity.FourierDensity;        // FourierDensity (ρ̃ for the §3 order-parameter diagnostic)
import qchem.CompositeCD;                         // tComposite_CD (the polarized density = one composite over Up+Down blocks, V1.37)
import qchem.ChargeDensity.Factory;
import qchem.ChargeDensity.SeedCD;              // PolarizedSeedCD (the raw spin-SAD seed, for the sublattice gate)               // IrrepCD_Factory/PolarizedCD_Factory (fixed-density probe)
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH, GTH_PP (the PP model, for the matrix-trace probe)
import qchem.Calculation;                        // qchem::Calculation, CalcOptions (finite reference)
import qchem.AtomCalculation;                    // AtomCalculation, AtomType, BasisSetAccuracy (Slater/High pseudo-atom ref)
import qchem.Types;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Complex_BS;
using qchem::BasisSet::Gaussian::BasisSetData;

using namespace qchem::tests::gpw;   // the harness: Materials cells, gates, recipes, probes (IntegrationTests/GPW/Harness.C)


// (PolarizedRunKeepsItsSpin -- the Mn sextet asking for Kerker, 217 s, 27% of the suite -- was DELETED
//  2026-09-15 (doc/TestSuitePlan.md phase 1).  Its claim is a MIXER property: a polarized seed must get a
//  per-channel ρ̃ mixer, never a single-map one.  It is now pinned with no SCF at all by
//  src/ChargeDensity/tests/KerkerMix.C -- PolarizedSeedComposesPerChannel, PolarizedStepMovesEachChannelAtAlpha
//  and the QCHEM_SPINBLIND_KERKER negative control -- in 41 ms.  What else it asserted is covered by
//  GPW_MnBox.Γ_M6_Smear_eqFinite above (same cell, the polarized energy) and GPW_Mn2Box.Γ_Becke_Shub_Pol_Smear_KeepsOrder
//  below (order SUSTAINED through an SCF).)

// ============================ MnO rocksalt AFM-II (SymmetryUpgradePlan §7 step 7) ============================
// The campaign RUN (RunMnO: the rhombohedral AFM-II cell, the recipe, the anneal schedule, the FM arm and the
// ordering comparison) is CLIapps/gpwprobe's `mno` sub-command since 2026-09-15 (doc/TestSuitePlan.md §8).  What
// stays here are the SEED-LEVEL gates: no SCF, cheap, and the mirror they pin is exact.

// THE RAW SEED of the MnO AFM-II cell: are the two Mn sublattices actually equal and opposite?
//
// One magnetic species on one Wyckoff site means the two Mn are related by the (1/2,1/2,1/2) translation,
// so ANY valid magnetic solution -- seed included -- must have m(Mn1) = -m(Mn2).  Nothing checked that.
// PW_Mn2Box.Γ_Pol_SeedStaggered covers a DIFFERENT cell (simple cubic, 2 Mn, no O, NEUTRAL
// targets) and asserts G-space quantities plus GetTotalSpin()==0; m_stag = 1/2(m1-m2), the campaign's order
// parameter, is BLIND to the imbalance -- it reads +0.366 for (+0.37,-0.37) and for (0,-0.73) alike.
// Measured 2026-08-11: the density one Fock build downstream of this seed has m1 = -0.0001, m2 = -0.73.
// This test asks whether the seed ITSELF is lopsided, i.e. whether the defect is upstream of the SCF.
TEST(GPW_MnO, Γ_Pol_SeedMirror)
{
    const double a=8.40;
    const Material mno=qchem::Materials::Get("MnO_AFM2");   // the rhombohedral AFM-II cell: Mn +m at 0 (the CORNER atom), Mn -m at 1/2, O at 1/4, 3/4
    const std::shared_ptr<UnitCell> cellp=mno.cell;
    UnitCell& cell=*cellp;
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The seed the run actually uses: PolarizedSeedCD over this cell, with the IonicSAD targets
    // (Mn2+ => 5 of the q7 valence, O2- => 8 of the q6).  Built directly -- no SCF, no Hamiltonian.
    qchem::BasisSet::PlaneWave::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    qchem::ChargeDensity::PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);
    const auto* up=seedCD.GetChannel(Spin::Up);
    const auto* dn=seedCD.GetChannel(Spin::Down);

    const rvec3_t off(0.7,0,0), rMn1(0,0,0), rMn2(a,a,a);   // A*(1/2,1/2,1/2) = a(1,1,1); 0.7 bohr = the d peak
    const double m1=(*up)(rMn1+off)-(*dn)(rMn1+off);
    const double m2=(*up)(rMn2+off)-(*dn)(rMn2+off);
    std::cout << "[seed probe] m1=" << m1 << " m2=" << m2 << " m1+m2=" << m1+m2
              << "  N_up=" << up->GetTotalCharge() << " N_dn=" << dn->GetTotalCharge() << std::endl;

    EXPECT_NEAR(up->GetTotalCharge(), dn->GetTotalCharge(), 1e-8) << "the AFM seed carries no NET moment";
    EXPECT_GT(std::abs(m1), 0.1) << "the + sublattice must actually be polarized";
    EXPECT_GT(std::abs(m2), 0.1) << "the - sublattice must actually be polarized";
    EXPECT_NEAR(m1+m2, 0.0, 0.02*std::abs(m1))
        << "the two sublattices are related by a lattice translation: m1 must equal -m2";

    // THE MIRROR RELATION EVERYWHERE, not just at the two nuclei.  The AFM seed satisfies
    // rho_up(r) = rho_dn(r + t) with t = A*(1/2,1/2,1/2) = (a,a,a) at EVERY r, by construction.  Probed
    // here at points that straddle the CELL BOUNDARY, because that is where the Becke XC mesh wraps its
    // points into the home cell (kpt = r - A*n0) and where a real-space evaluation built from per-atom
    // recentred radials WITHOUT periodic images would pick up the wrong atom -- which, the two Mn carrying
    // OPPOSITE spins, swaps the channels rather than merely losing density.  That failure is invisible to
    // the grid-charge check (rho_up+rho_dn is untouched; only rho_up-rho_dn flips) and it is exactly the
    // observed symptom: the corner atom's moment dies while the mesh integrates the total to 1e-5.
    struct P { const char* what; rvec3_t r; };
    const std::vector<P> probes = {
        {"inside, +x of Mn1",        rvec3_t( 0.7, 0.0, 0.0)},
        {"OUTSIDE the cell, -x",     rvec3_t(-0.7, 0.0, 0.0)},   // physically beside Mn1; wraps
        {"OUTSIDE the cell, -xyz",   rvec3_t(-0.5,-0.5,-0.5)},
        {"far tail, -2 bohr",        rvec3_t(-2.0, 0.0, 0.0)},
    };
    for (const auto& p : probes)
    {
        const rvec3_t rt = p.r + rvec3_t(a,a,a);                  // the mirror partner
        const double u=(*up)(p.r), d=(*dn)(rt);
        std::cout << "[mirror] " << p.what << ": rho_up(r)=" << u << "  rho_dn(r+t)=" << d
                  << "  diff=" << u-d << std::endl;
        EXPECT_NEAR(u, d, 1e-6*std::max(1.0,std::abs(u)))
            << "AFM mirror broken at " << p.what << ": the seed is not properly periodic there";
    }

    // THE BATCH OVERLOAD vs THE SINGLE-POINT ONE.  Everything above used operator()(rvec3_t).  The XC mesh
    // samples a MATRIX-FREE seed through the BATCHED operator()(rvec3vec_t) instead -- SinglesDensitySampler::RhoPol
    // takes its cSpinResolved_CD branch for exactly this density -- so the batch path is what the first Fock
    // build actually sees, and nothing has ever checked the two agree.  They must, pointwise.
    rvec3vec_t batch(2*probes.size());
    for (size_t i=0;i<probes.size();++i) { batch[i]=probes[i].r; batch[probes.size()+i]=probes[i].r+rvec3_t(a,a,a); }
    const rvec_t bu=(*up)(batch), bd=(*dn)(batch);
    ASSERT_EQ(bu.size(), batch.size());
    for (size_t i=0;i<batch.size();++i)
    {
        const double su=(*up)(batch[i]), sd=(*dn)(batch[i]);
        std::cout << "[batch] r=(" << batch[i].x << "," << batch[i].y << "," << batch[i].z << ")"
                  << " up: batch=" << bu[i] << " single=" << su << " d=" << bu[i]-su
                  << " | dn: batch=" << bd[i] << " single=" << sd << " d=" << bd[i]-sd << std::endl;
        EXPECT_NEAR(bu[i], su, 1e-6*std::max(1.0,std::abs(su))) << "UP batch != single at probe " << i;
        EXPECT_NEAR(bd[i], sd, 1e-6*std::max(1.0,std::abs(sd))) << "DN batch != single at probe " << i;
    }
}


// THE v_xc SUBLATTICE MIRROR ON THE XC MESH -- the SymmetryUpgradePlan "NEXT ACTION" probe (2026-08-11).
// By elimination (seed, Becke weights, Kinetic/Vloc/Vnl, Phi tables all exonerated) the first-Fock-build
// mirror break must live in v_xc -- yet a pointwise LSDA functional "cannot" be site-dependent.  This probe
// resolves the contradiction by testing what the Fock build ACTUALLY consumes: the channel rasters
// SinglesDensitySampler::RhoPol hands the spin-native XC term, at the mesh's own points.  The mesh stores its
// points WRAPPED into the home cell (kpt = r - A*n0, MakePeriodicBeckeMesh) -- so a valid seed must satisfy
// rho_up(p_g) = rho_dn(p_g + t) with t = A*(1/2,1/2,1/2) AT EVERY STORED POINT, and (v_xc being pointwise
// in the channel pair) v_xc^up(p_g) = v_xc^dn(p_g + t).  The two Mn blocks' grids are exact t-translates of
// each other (same radial x angular template), so every point's mirror partner is itself a mesh point --
// found here by hashing wrapped fractional coordinates, no interpolation anywhere.
//
// WHY THE SEED GATE ABOVE COULD PASS WHILE THIS FAILS: its probe points sit within ~2 bohr of a HOME atom,
// where SeedCD's real-space rho(r) = Sum_atoms rho_atom(|r-R|) -- a sum with NO LATTICE IMAGES -- is
// dominated by an atom it actually contains.  A WRAPPED mesh point near a cell face reads its density from
// an IMAGE of an atom (for the CORNER Mn, 7 of the 8 octants of its density hump belong to images), which
// an image-less sum simply does not have.  That is site-specific (the centre Mn2 has no near-shell wrapped
// points), a rigid translation changes it (MNO_SHIFT), the uniform raster reproduces it (its corner-region
// points read image density too), and under the AFM staggering the missing hump is the MAJORITY channel of
// exactly one sublattice -- every recorded symptom.
TEST(GPW_MnO, Γ_Becke_Pol_SeedVxcMirror)
{
    const double a=8.40;
    const Material mno=qchem::Materials::Get("MnO_AFM2");   // the rhombohedral AFM-II cell: Mn +m at 0 (the CORNER atom), Mn -m at 1/2, O at 1/4, 3/4
    const std::shared_ptr<UnitCell> cellp=mno.cell;
    UnitCell& cell=*cellp;
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The seed the run uses (identical to GPW_MnO.Γ_Pol_SeedMirror).
    qchem::BasisSet::PlaneWave::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    qchem::ChargeDensity::PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);

    // The run's Becke recipe, shrunk (nR=20, GL-11) -- wrapping is generic, it needs no production density.
    qcMesh::MeshParams mp=qcMesh::BeckeXCParams(20, -1.0, 11);
    auto mesh=std::make_shared<const qcMesh::Mesh>(cell.CreateIntegrationMesh(mp));
    const rvec3vec_t& P=mesh->Points();
    const size_t N=P.size();
    ASSERT_GT(N, 0u);

    // The channel rasters EXACTLY as the Fock build gets them (RhoPol's cSpinResolved_CD seed branch).
    auto engine=SinglesEngineOver({mesh, {}});
    const rvec_t ru=engine->RhoPol(&seedCD, Spin::Up);
    const rvec_t rd=engine->RhoPol(&seedCD, Spin::Down);

    // Mirror-partner lookup: hash each point's wrapped fractional coords; partner(g) = the mesh index at
    // wrapped(frac(p_g) + (1/2,1/2,1/2)).  Quantized key + 27-neighbour probe rides out wrap roundoff.
    const double q=1e-9;
    auto wrapf=[](rvec3_t f){ f.x-=floor(f.x); f.y-=floor(f.y); f.z-=floor(f.z); return f; };
    std::map<std::tuple<long long,long long,long long>,size_t> at;
    std::vector<rvec3_t> F(N);
    for (size_t g=0; g<N; g++)
    {
        F[g]=wrapf(cell.ToFractional(P[g]));
        at[{llround(F[g].x/q),llround(F[g].y/q),llround(F[g].z/q)}]=g;
    }
    auto partner=[&](size_t g)->long
    {
        const rvec3_t fm=wrapf(F[g]+rvec3_t(0.5,0.5,0.5));
        const long long kx=llround(fm.x/q), ky=llround(fm.y/q), kz=llround(fm.z/q);
        for (long long dx=-1; dx<=1; dx++) for (long long dy=-1; dy<=1; dy++) for (long long dz=-1; dz<=1; dz++)
            if (auto it=at.find({kx+dx,ky+dy,kz+dz}); it!=at.end())
            {
                rvec3_t d=wrapf(F[it->second]-fm+rvec3_t(0.5,0.5,0.5))-rvec3_t(0.5,0.5,0.5);  // min-image
                if (norm(d)<5e-9) return long(it->second);
            }
        return -1;
    };

    // v_xc per point from the channel pair -- the same functionals the spin-native XC term applies, through
    // the same two-channel face.
    qchem::Hamiltonian::SlaterExchange ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vc;
    auto vxc=[&](double u, double d, const Spin& s)->double
    {
        if (u+d<=1e-12) return 0.0;                       // VWN's r_s/log guard; symmetric, mirror-safe
        return ex.GetVxc(u, d, s) + vc.GetVxc(u, d, s);
    };

    // Sweep: rho and v_xc mirror defects across the whole raster; localize the worst offenders.
    size_t nOrphan=0, nBad=0;
    double maxRho=0, maxV=0, maxOrphanWRho=0;
    size_t argRho=0; long argRhoJ=-1;
    const rvec_t& W=mesh->Weights();
    for (size_t g=0; g<N; g++)
    {
        const long j=partner(g);
        // An ORPHAN is legal but must be an eps-tail point: the free Becke builder's keep decisions are
        // bit-different between translated blocks, so an eps-BORDERLINE point can be kept on one atom and
        // dropped on its partner (the same mechanism the imposed builder's orbit-consistency pass drops).
        // What the quadrature sees of it is w*rho -- bound THAT, not the count.
        if (j<0) { nOrphan++; maxOrphanWRho=std::max(maxOrphanWRho, std::abs(W[g])*std::max(ru[g],rd[g])); continue; }
        const double dRho=std::max(std::abs(ru[g]-rd[j]), std::abs(rd[g]-ru[j]));
        const double dV  =std::max(std::abs(vxc(ru[g],rd[g],Spin::Up  )-vxc(ru[j],rd[j],Spin::Down)),
                                   std::abs(vxc(ru[g],rd[g],Spin::Down)-vxc(ru[j],rd[j],Spin::Up  )));
        if (dRho>1e-6) nBad++;
        if (dRho>maxRho) { maxRho=dRho; argRho=g; argRhoJ=j; }
        if (dV  >maxV) maxV=dV;
    }
    std::cout << "[vxc mirror] N=" << N << " orphans=" << nOrphan << " (max w*rho " << maxOrphanWRho
              << ") bad(rho>1e-6)=" << nBad
              << "  max|rho_up(r)-rho_dn(r+t)|=" << maxRho
              << "  max|vxc_up(r)-vxc_dn(r+t)|=" << maxV << std::endl;

    // Localize the worst point: whose density is it -- a HOME atom's, or an IMAGE's the seed cannot see?
    if (argRhoJ>=0)
    {
        std::vector<rvec3_t> R; std::vector<std::string> nm={"Mn1","Mn2","O1","O2"};
        for (auto atom : cell) R.push_back(atom->itsR);
        auto nearest=[&](const rvec3_t& r)
        {
            double best=1e300; std::string who;
            for (size_t ia=0; ia<R.size(); ia++)
                for (int i=-1;i<=1;i++) for (int jj=-1;jj<=1;jj++) for (int k=-1;k<=1;k++)
                {
                    const double d=norm(r-(R[ia]+cell.ToCartesian(rvec3_t(i,jj,k))));
                    if (d<best) { best=d; who=nm[ia]+((i||jj||k)?" IMAGE":""); }
                }
            return std::make_pair(best,who);
        };
        const auto [dg,wg]=nearest(P[argRho]);
        const auto [dj,wj]=nearest(P[size_t(argRhoJ)]);
        std::cout << "[vxc mirror] worst point r=("<<P[argRho].x<<","<<P[argRho].y<<","<<P[argRho].z
                  << ") nearest "<<wg<<" d="<<dg<<"  rho_up(r)="<<ru[argRho]<<" rho_dn(r)="<<rd[argRho]<<"\n"
                  << "[vxc mirror] partner     r=("<<P[size_t(argRhoJ)].x<<","<<P[size_t(argRhoJ)].y<<","
                  << P[size_t(argRhoJ)].z<<") nearest "<<wj<<" d="<<dj
                  << "  rho_up(r+t)="<<ru[size_t(argRhoJ)]<<" rho_dn(r+t)="<<rd[size_t(argRhoJ)]<<std::endl;
    }

    EXPECT_LT(maxOrphanWRho, 1e-8) << "a partnerless mesh point must be an eps-tail point (weight*rho "
                                      "below the Becke builder's eps-converged-series contract)";
    EXPECT_LT(maxRho, 1e-8) << "the seed's channel rasters break the sublattice mirror ON THE XC MESH -- "
                               "this is the site-dependent v_xc defect (SymmetryUpgradePlan WHERE-WE-LEFT-OFF)";
    EXPECT_LT(maxV, 1e-6) << "v_xc^up(r) != v_xc^dn(r+t) on the mesh: the first Fock build is fed a "
                             "mirror-broken exchange-correlation potential";
}


// SHUBNIKOV S3, end to end at the seed level (doc/SymmetryUpgradePlan.md §7 step 7): a MAGNETICALLY
// imposed MnO basis must star-average the channel pair under the Shubnikov group -- keeping the AFM
// staggering EXACTLY mirror-symmetric -- where the grey average would erase it.  Chain under test:
// MagneticDecoration (the seed's own species rule) -> GPWParams::siteSpins -> the factory's Shubnikov
// resolution -> CreateXCQuadrature (site-adapted invariant mesh + fold + sigma tags + flip-fixed zero
// flags) -> SinglesDensitySampler::RhoPol's (rho,m) split.
TEST(GPW_MnO, Γ_Shub_Pol_SeedDecoration)
{
    namespace L3=BasisSet::Lattice;
    using namespace qchem::ChargeDensity;
    const double a=8.40;
    const Material mno=qchem::Materials::Get("MnO_AFM2");   // the rhombohedral AFM-II cell: Mn +m at 0 (the CORNER atom), Mn -m at 1/2, O at 1/4, 3/4
    const std::shared_ptr<UnitCell> cellp=mno.cell;
    UnitCell& cell=*cellp;
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The decoration, by the seed's own resolution: Mn2+ (d5 pair) = +/-1, O2- (closed shell) = 0.
    auto st=lat.GetStructure();
    std::vector<int> spins=MagneticDecoration(st.get(), "LDA", IonicSADTargets(st.get(), "LDA"));
    ASSERT_EQ(spins.size(), 4u);
    EXPECT_EQ(spins[0], +1); EXPECT_EQ(spins[1], -1);
    EXPECT_EQ(spins[2],  0); EXPECT_EQ(spins[3],  0);

    // The magnetically IMPOSED GPW basis (coarse everything: this gate tests symmetry, not accuracy).
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR),
        L3::GPWParams{.densityEcut=8.0, .imposeSymmetry=true, .siteSpins=spins}));
    BasisSet::FitQuadrature q = bs->CreateXCQuadrature(st.get(), qcMesh::BeckeXCParams(20, -1.0, 11));
    ASSERT_EQ(q.NumSpinOps(), 24u) << "the Shubnikov group of the AFM-II cell has 24 ops (12+12)";
    size_t nFlip=0; for (auto s : q.GetSpinOps()) if (s==Symmetry::SpinAction::Flip) nFlip++;
    EXPECT_EQ(nFlip, 12u);
    ASSERT_EQ(q.NumFlipFixed(), q.GetMesh()->size());   // (FoldedMesh's ctor now checks this too)

    // The seed (the same PolarizedSeedCD the run uses), through the engine's channel-pair projector.
    qchem::BasisSet::PlaneWave::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);

    auto engine=SinglesEngineOver(q);
    const rvec_t up=engine->RhoPol(&seedCD, Spin::Up);
    const rvec_t dn=engine->RhoPol(&seedCD, Spin::Down);

    // (a) The staggering SURVIVES the imposed projector: the magnetization raster keeps its full scale.
    double maxM=0; for (size_t g=0; g<up.size(); g++) maxM=std::max(maxM, std::abs(up[g]-dn[g]));
    EXPECT_GT(maxM, 0.1) << "the SHUBNIKOV average must PRESERVE the AFM staggering";

    // (b) ...and is EXACTLY mirror-antisymmetric on the raster: m must vanish at every flip-fixed point,
    // and the weighted net moment must be zero to machine precision (the signed projector's guarantees).
    const rvec_t& W=q.GetMesh()->Weights();
    double net=0, worstFixed=0;
    for (size_t g=0; g<up.size(); g++)
    {
        net += W[g]*(up[g]-dn[g]);
        if (q.IsFlipFixed(g)) worstFixed=std::max(worstFixed, std::abs(up[g]-dn[g]));
    }
    // The projector pairs every orbit's members +/- exactly; the residual is the ULP-level asymmetry of
    // the site-adapted builder's partner WEIGHTS (op-image copies, equal only to roundoff -- measured
    // 1.7e-9 over 6000 points), not a projector defect.
    EXPECT_LT(std::abs(net), 1e-7) << "an imposed AFM pair carries ZERO net moment by construction";
    EXPECT_LT(worstFixed, 1e-12) << "m must vanish exactly at the flip-fixed mesh points";

    // (c) The GREY control: the SAME mesh and fold with the sigma tags withheld = the historical
    // per-channel spatial average, which maps +m sites onto -m sites and ERASES the order.  This is the
    // unit-level half of the S4 negative control ("imposing the grey group kills m_stag").
    auto grey=SinglesEngineOver(qcMesh::FoldedMesh(q.GetMesh(), q.GetFold()));
    const rvec_t gup=grey->RhoPol(&seedCD, Spin::Up);
    const rvec_t gdn=grey->RhoPol(&seedCD, Spin::Down);
    double maxGrey=0; for (size_t g=0; g<gup.size(); g++) maxGrey=std::max(maxGrey, std::abs(gup[g]-gdn[g]));
    EXPECT_LT(maxGrey, 1e-10) << "the grey average must ERASE the staggering -- if it does not, the "
                                 "Shubnikov machinery is not actually load-bearing";

    // (d) The TOTAL density is identical through both engines (the even channel is sigma-blind).
    double dTot=0; for (size_t g=0; g<up.size(); g++) dTot=std::max(dTot, std::abs((up[g]+dn[g])-(gup[g]+gdn[g])));
    EXPECT_LT(dTot, 1e-10) << "sigma must not touch the total density";
}
