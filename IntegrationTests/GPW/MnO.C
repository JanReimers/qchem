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
//   GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Long      (LONG: ~2.5 min; `ctest -L long` -- the converged AFM-II gate)
//   GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long      (LONG: two ~6 min arms; the shell-averaged DFT+U oracle gate, VA span)

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
#include <filesystem>    // the deck run's scratch directory
#include <nlohmann/json.hpp>  // the deck the AFM arm is run from

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
import qchem.RunPolicy;                         // SetRunPolicy / SolidCalcOptions::policy -- the declared-deviation A/B hatch (N5)
import qchem.SolidCalculation;                    // the NAMED periodic facade (Step 4 3/3)
import qchem.Deck;                            // deck::Run -- the AFM arm below is run FROM A DECK (D-ENV step 6d.3)
import qchem.Tests.GPW_Harness;                   // THE HARNESS (IntegrationTests/GPW/Harness.C): Materials cells, gates, recipes, the XC probes
import qchem.Materials;                           // Materials::Get -- the cells come from src/Structure/Data/materials.json (row MD)
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

// ============================ THE CONVERGED AFM-II GATE (back as a cell, 2026-09-15) ============================
// The old DISABLED_MnO_AFM2_RhombohedralGamma carried two things: the CAMPAIGN INSTRUMENT (a cell with three
// geometry discriminators, a knob-driven recipe, the FM arm, the ordering comparison -- now `gpwprobe mno`) and
// the GATE hidden inside it: does the converged AFM-II state reproduce?  This is the gate, and nothing else.
// THE RECIPE (re-judged 2026-09-20 under doc/Benchmark.md rule 3f -- the DECK-SHAPED loop): IonicSAD seed
// (Mn2+ d^5 + O2-), pivoted Cholesky at 1e-4 (cond(S)~7e8), Kerker G0=1 against the low-G charge-transfer
// slosh with ONE density-side history (Pulay depth 8 after 5 priming steps) and NOTHING on the Fock side, no
// MOM, kT=5e-3 riding the open d manifold, alpha=0.45, CP2K's convergence measure max|dD| < 1e-6 -- and the
// SHUBNIKOV group of the seed's decoration imposed (S3), which is what holds the staggering exactly mirrored
// through the loop.  The 2026-09-15 recipe (Fock DIIS->GDM Ladder + delayed MOM at 1e-5 on the Kerker
// residual) "converged" in 42 iterations while its density matrix was still moving at 3e-3 per iteration:
// two histories extrapolating each other's output, and MOM rotating a smeared degenerate d frontier (the
// +U gate below found it; OpenWork step 5).  The machinery tokens are elided from the name as NaF's are: they
// are this material's production recipe, not what the test is about.
// MEASURED 2026-09-20 (this recipe): 17 iterations to max|dD| 6e-7, 2m33s wall, Etot=-61.41454675, integrated
// site moment 4.451 e -- the same energy to 2e-7 Ha as the old recipe's 42 iterations / 5m20s / -61.41454697
// (pin 10: an anchor is re-judged against an INDEPENDENT ROUTE, and a different loop reaching the same number is
// one).  |m-tilde(q_AFM)| Omega/2 = 3.13 e.  CP2K's AFM-II oracle is
// -61.470570 (deck IntegrationTests/CP2K/mno_afm2_gpw_sr.inp): the 56 mHa gap is the banked ordering/d-selective
// offset (doc/SymmetryUpgradePlan.md §7, doc/SphericalLatticePlan.md), so this pins OUR number as a did-E-move
// anchor and states the oracle beside it.  LONG (ruling 5: over the 60 s budget), hence the `_Long` suffix: it runs in the
// `ctest -L long` tier (excluded by `ctest -LE long`), and in the TestMate tree like any other test.
TEST(GPW_MnO, Γ_Shub_Pol_Smear_Anchor_Long)
{
    const Material mno=qchem::Materials::Get("MnO_AFM2");   // Mn +m at 0, Mn -m at 1/2 (the decoration), O at 1/4, 3/4
    const Lattice_3D lat=LatticeOf(mno);
    SolidCalcOptions o=OptionsFor(mno, "MnO AFM-II Gamma (imposed)");
    o.multiplicity=1;                                       // the explicit two-channel singlet: nUp=nDn=13
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.accelerator=qchem::SCFAccelerators::Type::Null;      // the deck-shaped loop: the history is on the density side
    o.imposeSymmetry=true;                                  // S3: the Shubnikov group of the declared ordering
    SCFParams par=Gates(200, 1e-6, 1e30);
    par.Δρmeasure=SCFParams::Measure::MaxΔD;               // CP2K's EPS_SCF measure (rule 3f)
    par.PulayDepth=8; par.PulayStart=5;
    par.StartingRelaxRo=0.45; par.KerkerG0=1.0;
    par.UseMOM=false;
    par.SmearingkT=5e-3;
    Trace trace; o.onIteration=trace.Observer();
    GpwReport report("MnO "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*mno.cell, BasisSetData::VALENCE_LOWQ_SR), o, par);
    trace.Print(o.label, /*polarized*/true);
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 26.0, 1e-6);
    EXPECT_NEAR(R->Energy(), -61.41455, 2e-3) << "did-E-move anchor (CP2K AFM-II oracle -61.470570: the banked 56 mHa ordering offset)";
    ASSERT_FALSE(trace.rows.empty());
    EXPECT_GT(std::abs(trace.rows.back().order), 4.0) << "the INTEGRATED Mn site moment must be d^5-scale at convergence (measured 4.45)";
    EXPECT_NE(R->SpinDensity(), nullptr) << "a polarized run must hand back m(r)";
}

// ============================ THE SHELL-AVERAGED DFT+U ORACLE GATE (programme step 5.1) ============================
// CP2K's +U oracle (doc/Records/CP2Kresults.md, 2026-09-19): the banked VA deck + `&DFT_PLUS_U L 2,
// U_MINUS_J [eV] 4.0` on both Mn kinds, PLUS_U_METHOD LOWDIN -> Etot -60.68597088 Ha (was -61.30332518 without
// U: a +0.61735 Ha shift), DFT+U energy 0.60950923 Ha, Mulliken moments +-4.765.  What is compared is the
// U-INDUCED SHIFT and E_U itself, not the absolute energy: the two codes already sit 100 mHa apart on this span
// (doc/Benchmark.md section 5c, the XC-grid/fit deviations), and that offset is U-independent to first order.
// The manifold is CP2K's definition, read off src/dft_plus_u.F: EVERY d shell on the Mn -- all 8 contractions,
// a 40x40 Loewdin block per Mn per spin -- which qchem's Hubbard_U reproduces (Columns).  AND THE FORM IS CP2K's
// (found 2026-09-20 when the Dudarev form gave E_U=0.080 against the oracle's 0.610): dft_plus_u.F keeps ONLY
// THE DIAGONAL of the block (`IF (isgf == jsgf)`), i.e. the 40 Loewdin POPULATIONS, no eigen-decomposition --
// QCHEM_U_EIGEN=0 here, the CP2K_COMPAT member the user ruled it into (2026-09-20) -- so the two E_U's are the
// same functional of two densities.  The Dudarev form on
// the same recipe, MEASURED 2026-09-20 (self-consistent, U=4 eV): E_U=0.0801 Ha, dE=+0.1034 Ha, occupation
// eigenvalues near-integer (N_maj=4.80, max lambda=0.997, N_min=0.25) -- the physical number, 7.6x below the
// population form's; it is banked in doc/Records/CP2Kresults.md, not asserted here.
// That is why this needs the VA span (the 8-zeta d file CP2K holds
// function-for-function) AND the spherical view: a Cartesian d has six components with the s-contaminant among
// them, and Hubbard_U refuses it (the Loewdin block would be 48x48 with a contaminant in the manifold).
// The RECIPE is the Shub anchor's (above), which is the production one for this cell.  Two arms, ~6 min each
// (the VA span is the 6m38s benchmark row), hence `_Long` like the anchor (`ctest -L long`).
// QCHEM_U_TRACE=1 prints the per-refresh occupation line (N per manifold, max lambda, E_U) to compare against
// CP2K's own &PRINT &PLUS_U table in bench_MnO_AFM2_VA_plusU_cp2k.log.
TEST(GPW_MnO, Γ_U_Shub_Pol_Smear_CP2K_Long)
{
    const Material mno=qchem::Materials::Get("MnO_AFM2");   // sites 0,1 = Mn (+m, -m); 2,3 = O
    const Lattice_3D lat=LatticeOf(mno);
    {   // the manifold indices are a claim about the material row: check it, do not assume it
        int i=0; mno.cell->ForEachSite([&](int Z, const rvec3_t&, bool){ if (i<2) EXPECT_EQ(Z,25) << "site "<<i<<" must be Mn"; i++; });
    }
    auto basis=[&]{ return BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(
        std::shared_ptr<const Real_BS>(BasisSet::Gaussian::Factory(BasisSetData::VALENCE_LOWQ_VA, mno.cell.get(),
                                       BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian))); };
    // Each arm's SolidCalculation must OUTLIVE its Result (the result points into the calculation), so the
    // arms are held here and the lambda only builds them.
    std::vector<std::unique_ptr<qchem::SolidCalculation>> arms;
    auto run=[&](double U_eV, const std::string& label)
    {
        SolidCalcOptions o=OptionsFor(mno, label);
        o.multiplicity=1;
        o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
        o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
        // THE CP2K-SHAPED LOOP (user, 2026-09-20: "we need 1) very close converged energies 2) similar SCF
        // convergence rates 3) similar runtime and RAM").  The deck: BROYDEN_MIXING (Kerker BETA 1.5, NBUFFER 8)
        // on rho-tilde, diagonalise, EPS_SCF 1e-6 on max|dP|, MAX_SCF 200, no Fock-side extrapolation, no MOM.
        // Ours, matched: ONE density-side history (Kerker-preconditioned Pulay, depth 8) and NOTHING on the Fock
        // side -- the anchor's recipe (Ladder DIIS + MOM) on top of that history left max|dD| at 3e-3 for 200
        // iterations while the G-space residual read "converged" at 43: Fock-DIIS and density-Pulay extrapolate
        // each other's output, and MOM hands a smeared degenerate d frontier unequal occupations that swap every
        // iteration (D rotates at constant rho).  With the deck's loop shape and its measure: 22 iterations to
        // max|dD| 6e-8 against the deck's 44 (probe, `gpwprobe mno`, MNO_ACC=Null MNO_MOM=0 MNO_PULAY=8).
        o.accelerator=qchem::SCFAccelerators::Type::Null;
        o.imposeSymmetry=true;
        o.policy.hubbardEigen=false;   // CP2K's form, a declared deviation: stated in the options (D-ENV 6b), not in the environment
        if (U_eV>0.0) o.hubbard={HubbardU(0,2,U_eV), HubbardU(1,2,U_eV)};
        SCFParams par=Gates(200, 1e-6, 1e30);
        par.Δρmeasure=SCFParams::Measure::MaxΔD;             // CP2K's EPS_SCF measure (doc/Benchmark.md rule 3f)
        par.PulayDepth=8; par.PulayStart=5;                   // the density history, as NBUFFER 8
        par.StartingRelaxRo=0.45; par.KerkerG0=1.0;
        par.UseMOM=false;
        par.SmearingkT=5e-3;
        Trace trace; o.onIteration=trace.Observer();
        GpwReport report("MnO "+o.label, par.Verbose);
        arms.push_back(std::make_unique<qchem::SolidCalculation>(lat, basis(), o, par));
        trace.Print(o.label, /*polarized*/true);
        return arms.back()->Result();
    };
    // THE FORM is a declared deviation, so it is set the way the N5 hatch sets one: for this test only.
    // (THE FORM -- CP2K's diagonal populations -- is stated in the run's options below: o.policy.hubbardEigen=false.)
    auto r0=run(0.0, "MnO AFM-II VA sph Gamma (imposed)");
    auto rU=run(4.0, "MnO AFM-II VA sph Gamma (imposed, U=4 eV on Mn d)");
    ASSERT_TRUE(r0) << Why(r0);
    ASSERT_TRUE(rU) << Why(rU);
    EXPECT_NEAR(r0->TotalCharge(), 26.0, 1e-6);
    EXPECT_NEAR(rU->TotalCharge(), 26.0, 1e-6);
    const double E0=r0->Energy(), EU=rU->Energy(), dE=EU-E0, E_U=rU->EnergyTerms()["E_U"];
    const double cp2k_dE=-60.68597087953861-(-61.30332518), cp2k_EU=0.60950923227552;
    std::cout<<std::setprecision(10)<<"[MnO +U] E(0)="<<E0<<"  E(U)="<<EU<<"  dE="<<dE<<" (CP2K "<<cp2k_dE
             <<")  E_U="<<E_U<<" (CP2K "<<cp2k_EU<<")  dE-E_U="<<dE-E_U<<" (CP2K "<<cp2k_dE-cp2k_EU<<")"<<std::endl;
    EXPECT_GT(E_U, 0.0);
    EXPECT_NEAR(E_U, cp2k_EU, 0.05) << "the same population functional (all 8 d shells, diagonal) of two densities that agree to 0.1 Ha";
    EXPECT_NEAR(dE,  cp2k_dE, 0.05) << "the U-induced shift; the absolute offset between the codes is U-independent";
}


// ============================ THE ORDERING GATE (was `gpwprobe mno`'s last check; D-ENV step 6d.3) ============================
// "AFM-II is the LSDA ground-state ordering of MnO": E(AFM-II) < E(FM), the same cell, the same recipe.  It used to be a PASS/FAIL line in a
// command-line probe driven by environment variables; the user ruled (2026-10-04) that an important gate is a hard-coded integration test.
// The AFM arm is run FROM A DECK (deck::Run -- which also gates the deck system on the production recipe: the anchor above pins the SAME
// energy built by hand, to its tolerance); the FM arm is built here, in the test, because FM-vs-AFM is a property of THIS test and not of the
// calculation framework: the same Bravais cell with no spin flip, multiplicity 11 (two d^5 Mn).  LONG: two arms of ~7 min.
// ⚠ DISABLED_ because the CLAIM IS CURRENTLY FALSE (CLAUDE.md: a real claim currently failing is an open tracker row, not a green test):
// MEASURED 2026-10-04, E_AFM=-61.414547 vs E_FM=-61.452697 -- the FM arm is 38 mHa BELOW (doc/OpenWork.md §3, "MnO ordering at Γ/SR").  The
// deck-vs-anchor assertion and the FM charge assertion below passed; only the ordering EXPECT failed.  Re-enable when the ordering is right.
TEST(GPW_MnO, DISABLED_Γ_Shub_Pol_Smear_Ordering_Long)
{
    namespace fs=std::filesystem;
    const nlohmann::json deck={
        {"structure","MnO_AFM2"}, {"basis",{{"data","VALENCE_LOWQ_SR"}}},
        {"solid",{{"multiplicity",1},{"seed","IonicSAD"},{"ortho","CholeskyPivoted"},{"orthoTol",1e-4},{"accelerator","Null"},{"imposeSymmetry",true}}},
        {"scf",{{"NMaxIter",200},{"minDeltaRho",1e-6},{"deltaRhoMeasure","MaxDeltaD"},{"minDeltaE",1e30},{"minDeltaFD",1e30},{"minVirial",1e30},{"minFD",1e30},
                {"startingRelaxRo",0.45},{"mergeTol",1e-4},{"pulayDepth",8},{"pulayStart",5},{"kerkerG0",1.0},{"useMOM",false},{"smearingkT",5e-3}}}};
    qchem::deck::RunSpec spec; qchem::deck::FromJson(deck, spec);
    const fs::path dir=fs::temp_directory_path()/("qchem_it_mno_ordering_"+std::to_string(::getpid()));
    fs::remove_all(dir);
    qchem::deck::Provenance pv; pv.codeVersion="it";
    const auto afm=qchem::deck::Run(spec, pv, dir);
    ASSERT_TRUE(afm.converged) << afm.summary;
    ASSERT_TRUE(afm.energy.has_value());
    EXPECT_NEAR(*afm.energy, -61.41455, 2e-3) << "the deck run reproduces the hand-built anchor (GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Long)";

    // the FM arm: same cell, NO spin flip, multiplicity 11, the same recipe
    auto cellp=std::make_shared<UnitCell>(BravaisCell(Bravais::CubicF, {.a=8.4}, Matrix3D<int>(0,1,1, 1,0,1, 1,1,0)));
    cellp->AddAtom(25, {0.0,0.0,0.0}, false);  cellp->AddAtom(25, {0.5,0.5,0.5}, false);
    cellp->AddAtom(8,  {0.25,0.25,0.25});     cellp->AddAtom(8,  {0.75,0.75,0.75});
    Lattice_3D lat(*cellp, ivec3_t(1,1,1));
    qchem::deck::RunSpec resolved=spec; qchem::deck::Resolve(resolved);     // fills Nelec / species from the material, as the AFM arm's run did
    SolidCalcOptions o=resolved.solid; o.multiplicity=11; o.label="MnO FM Gamma (imposed)";
    const SCFParams par=spec.scf;
    GpwReport report("MnO FM", false);
    qchem::SolidCalculation fm(lat, MakeBasisLowQ(*cellp, BasisSetData::VALENCE_LOWQ_SR), o, par);
    auto R=fm.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 26.0, 1e-6);
    std::cout << "[MnO ordering] E_AFM=" << std::setprecision(10) << *afm.energy << "  E_FM=" << R->Energy()
              << "  dE=" << (R->Energy()-*afm.energy)*1000 << " mHa" << std::endl;
    EXPECT_LT(*afm.energy, R->Energy()) << "AFM-II is the LSDA ground-state ordering";
    fs::remove_all(dir);
}
