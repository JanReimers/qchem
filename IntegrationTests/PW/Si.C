// File: IntegrationTests/PW/Si.C  Si diamond on the plane-wave basis (the framework anchors + the prototype-loop reference).
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(PW_<Material>, <k>_[tokens]_<Claim>); `Prototype` marks the standalone
// loop.  `scripts/testgrid` renders the coverage table from these names.
//
//   PW_Si.Γ_Prototype_Converges
//   PW_Si.Γ_eqPrototype
//   PW_Si.Γ_Anchor
//   PW_Si.k222_Prototype_Converges
//   PW_Si.k222_Anchor

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
#include "gtest/gtest.h"

import qchem.Tests.PW_Harness;   // the PW drivers + the direct-grid oracles + PWFixture
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;
import qchem.BasisSet.PlaneWave.Evaluators;   // PW_Grid_Evaluator -- a UNIT TEST may reach the internal
import qchem.BasisSet.Lattice.BasisSet;   // Factory(Type::PW, lat, Ecut, loc, nl) -> Complex_BS*
import qchem.ScalarFunction;                 // ScalarFunction<double> -- arg of the moved-here field oracles
import qchem.Mesh;                           // qcMesh::MeshParams (Vee_Hartree's fit-basis factory arg; ignored)
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Ewald;          // EwaldEnergy (ion-ion Madelung term -> physical total energy)
import qchem.Types;          // dcmplx, ivec3_t, rvec_t, mat_t, chmat_t
import qchem.Blaze;          // mat_t<dcmplx>
import qchem.Math;           // Pi
import qchem.Hamiltonian.Internal.ExFunctional;     // the validated LDA functional interface
import qchem.Hamiltonian.Internal.SlaterExchange;   // Dirac exchange (alpha=2/3), eps_x = 3/4 v_x
import qchem.Hamiltonian.Internal.VWN_Correlation;  // VWN5 correlation (validated vs libxc)
import qchem.Hamiltonian.Internal.PWTerms;          // Ven_PP_Short/Long, Vee_Hartree (dcmplx Hamiltonian terms)
import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
import qchem.Hamiltonian;                           // cStatic_HT / cDynamic_HT aliases (public term interfaces)
import qchem.Hamiltonian.Internal.Hamiltonian;      // cHamiltonianImp (the dcmplx Hamiltonian = sum of terms)
import qchem.Hamiltonian.Internal.Hamiltonians;     // Ham_PW_DFT (the assembled plane-wave LDA KS Hamiltonian)
import qchem.Hamiltonian.Internal.IonIon;           // IonIon<double> (ion-ion term: pair sum / Ewald via isFinite)
import qchem.Hamiltonian.Internal.Kinetic;          // Kinetic<dcmplx> (the shared kinetic term)
import qchem.Pseudopotential.GTH_Potentials;    // GetGTH (CP2K GTH/HGH database reader)
import qchem.Energy;                                // EnergyBreakdown
import qchem.ChargeDensity.Imp.IrrepCD;             // PeriodicIrrepCD<dcmplx> (the periodic leaf -- 3c-2b)
import qchem.ChargeDensity.SeedCD;                  // SeedCD + PolarizedSeedCD (the spin-SAD staggering gate)
import qchem.Fitting.FunctionFitter;                // Factory / ProjectedDensity_G / FunctionFitter_Density (item B)
import qchem.BasisSet.G_FieldEvaluator;             // the grid-engine seam (GridPoints/RhoOnGrid/Integral) for the item-K probe
import qchem.BasisSet.GMap;                             // ΔG_Map (the G-space coefficient map)
import qchem.BasisSet.Orbital_DFT_IBS;                      // cFIT_CD_ABS (the ortho density-fit basis face)
import qchem.Symmetry.Irrep;                        // Irrep
import qchem.LASolver;                              // complex Hermitian eigensolver
import qchem.Structure;                             // Molecule, Atom (the Si diamond basis)
import qchem.Matrix3D;                              // Matrix3D<double> (the FCC cell matrix)
import qchem.SCFIterator;                           // cSCFIterator (the real framework SCF driver)
import qchem.SCFParams;                             // SCFParams
import qchem.WaveFunction;                          // cWaveFunction (read view of the converged state)
import qchem.ElectronConfiguration;                 // ElectronConfiguration base
import qchem.ElectronConfiguration.Crystal;         // Crystal_EC (single-k Bloch configuration)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // SCFAcceleratorNull (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;      // SCFAcceleratorDIIS (scalar-agnostic manager)
import qchem.BasisSet.Internal.BasisSetImp;         // BasisSetImp<dcmplx> (single-block BasisSet container)

using namespace qchem;
using namespace qchem::tests::pw;

using BasisSet::PlaneWave::PlaneWave_IBS;
using BasisSet::PlaneWave::PW_Grid_Evaluator;   // internal grid evaluator (unit test may reach it directly)
using Pseudopotential::HGH_LocalPotential;
using Pseudopotential::HGH_SeparablePotential;
using Pseudopotential::GetGTH;
using Pseudopotential::GTH_PP;

namespace
{

// Real material: silicon, diamond structure, GTH-LDA q4 pseudopotential (local + KB nonlocal).
// FCC PRIMITIVE cell (2 Si, 8 valence electrons) at Gamma.  G=0 is dropped (neutralising background)
// so the ABSOLUTE total energy is shifted and not meaningful here (per the plan) -- we validate the
// SELF-CONSISTENT loop instead: it converges, the charge sum rule gives exactly 8 valence electrons,
// and the band-sum and direct total energies agree at the fixed point.  A modest Ecut keeps the direct
// (non-FFT) transforms fast; the grid is auto-sized to resolve the basis difference set.
TEST(PW_Si, Γ_Prototype_Converges)
{
    const double a=10.26, h=0.5*a;                         // Si lattice constant ~5.43 A in bohr
    Matrix3D<double> A(0.0,h,h,  h,0.0,h,  h,h,0.0);       // FCC primitive: cols = a/2 (011),(101),(110)
    UnitCell          cell(A);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    const double      Omega=cell.GetCellVolume();          // = a^3/4
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    Molecule si;                                           // 2-atom diamond basis
    si.Insert(new Atom(14, rvec3_t(0,0,0)));
    si.Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));

    GTH_PP                 siPP=GetGTH("Si","LDA",4);          // CP2K GTH-LDA q4, from the database
    const HGH_LocalPotential&     loc=siPP.local;
    const HGH_SeparablePotential& nl =siPP.nonlocal;
    chmat_t Vext = pw.MakeSpeciesFieldMatrix(&si, loc, qchem::BasisSet::FieldRange::Full) + pw.MakeProjectorMatrix(&si,nl);

    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};

    const int Nval=8;                                      // 2 Si x Zion 4
    int       m=MaxGComponent(pw);
    ivec3_t   Ng(4*m+1, 4*m+1, 4*m+1);                     // resolves the difference set (no aliasing)

    SCFResult R=RunSCF(pw, recip.GetCell(), Omega, Vext, Nval, Ng, vxcOf, epsOf, 0.4, 1e-8, 400);

    // Ion-ion (Ewald) energy of the two Si4+ cores: turns the electronic energy (which drops the G=0
    // potential and has no ion-ion term) into a physical, NEGATIVE total.  The Ewald cell carries the
    // same FCC geometry and the diamond basis given in fractional coordinates; charges = Zion = 4.
    UnitCell ecell(A);
    ecell.AddAtom(14, rvec3_t(0.0,0.0,0.0));
    ecell.AddAtom(14, rvec3_t(0.25,0.25,0.25));
    double Eii  = EwaldEnergy(ecell, rvec_t{4.0,4.0});

    // G=0 alignment of the local pseudopotential: the dropped G=0 potential carries a finite
    // electron-ion shift E_alpha = (N/Omega) Sum_a alpha_a, alpha_a = integral[V_loc^a + Zion/r].
    // (Two Si atoms, same species.)  This is the last G=0 piece needed for the absolute total energy.
    double Ealpha = (double(Nval)/Omega) * 2.0 * loc.FormFactorG0(14);
    double Etot   = R.Etot_direct + Eii + Ealpha;

    std::cout << "[Si] nG="<<pw.GetNumFunctions()<<" grid="<<(4*m+1)<<"^3 iters="<<R.iters
              << " converged="<<R.converged << "\n  Ekin="<<R.Ekin<<" Eext="<<R.Eext
              << " E_H="<<R.EH<<" E_xc="<<R.Exc << "\n  Etot(elec)="<<R.Etot_direct
              << " E_ion-ion(Ewald)="<<Eii << " E_alpha(G=0)="<<Ealpha
              << "\n  Etot(physical)="<<Etot << std::endl;

    ASSERT_TRUE(R.converged);
    EXPECT_NEAR(Omega*std::real(RhoAt(R.rho,ivec3_t(0,0,0))), double(Nval), 1e-6);  // 8 valence e-
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-5);                                  // stationary fixed point
    EXPECT_LT(Etot, 0.0);   // ion-ion Madelung + G=0 alignment make the total energy negative
}


// THE STAGE-3 PAYOFF: silicon (Gamma) self-consistent DFT run entirely through the FRAMEWORK objects --
// a cHamiltonianImp summing the PW Kohn-Sham terms (kinetic + external PP + Hartree + Dirac + VWN),
// an PeriodicIrrepCD<dcmplx> density built from the complex orbitals, and the framework energy bookkeeping --
// reproducing the standalone prototype's Si-Gamma result (Etot=1.468, 8 valence electrons).  A thin SCF
// driver stands in for the (not-yet-complexified) WaveFunction/SCFIterator orchestration.
TEST(PW_Si, Γ_eqPrototype)
{
    using namespace qchem::Hamiltonian;
    const double a=10.26, h=0.5*a;
    Matrix3D<double> Amat(0.0,h,h,  h,0.0,h,  h,h,0.0);
    UnitCell          cell(Amat);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    auto si=std::make_shared<Molecule>();
    si->Insert(new Atom(14, rvec3_t(0,0,0)));
    si->Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));
    GTH_PP                 siPP=GetGTH("Si","LDA",4);          // CP2K GTH-LDA q4, from the database
    const HGH_LocalPotential&     loc=siPP.local;
    const HGH_SeparablePotential& nl =siPP.nonlocal;

    // Framework Hamiltonian: a cHamiltonianImp summing the PW Kohn-Sham terms.  The external term
    // owns the pseudopotential model (the pseudo-wall) and assembles it through the basis.
    cHamiltonianImp ham(SpinGroup::UnPolarized);
    ham.Add(new Kinetic<dcmplx>);
    ham.Add(new Ven_PP_Short(si, &loc));                                // electron-ion SHORT-range local
    ham.Add(new Ven_PP_NonLocal(si, &nl));                              // KB separable projectors
    // The LONG-range core-charge V_long is its OWN term (the CP2K local-PP split, doc/GPWPlan.md 0e-PP);
    // a manual Hamiltonian must add it explicitly or it drops V_long entirely.
    ham.Add(new Ven_PP_Long(si, &loc));                                 // electron-ion LONG-range
    ham.Add(new Vee_Hartree(Vee_Hartree::fbs_t(pw.CreateCDFitBasisSet(nullptr, qcMesh::MeshParams{}))));
    ham.Add(NewPWXC(pw, std::make_shared<SlaterExchange> (2.0/3.0)));   // Dirac exchange
    ham.Add(NewPWXC(pw, std::make_shared<VWN_Correlation>()));          // VWN5 correlation

    const int Nelec=8, Nocc=Nelec/2;
    size_t  n=pw.GetNumFunctions();
    Irrep   irr=pw.GetIrrep(Spin::None);
    chmat_t S=pw.MakeOverlap();

    // Seed a UNIFORM density (Hartree+XC present from iteration 0, as real PW codes do): D = (N/n) I.
    hmat_t<dcmplx> Dprev=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++) Dprev(i,i)=double(Nelec)/double(n);
    auto cd=std::make_unique<qchem::ChargeDensity::PeriodicIrrepCD<dcmplx>>(Dprev, &pw, irr);

    bool converged=false; double Eprev=1e30;
    qchem::EnergyBreakdown E;
    for (int it=0; it<400; it++)
    {
        chmat_t H=ham.GetMatrix(&pw, Spin::None, cd.get());      // sums the framework terms for the current density
        LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
        las->SetBasisOverlap(S);
        auto sol=las->Solve(H);
        delete las;
        mat_t<dcmplx> U=std::get<0>(sol); rvec_t e=std::get<1>(sol);

        std::vector<size_t> idx(e.size());                       // occupy the Nocc lowest bands
        for (size_t i=0;i<idx.size();i++) idx[i]=i;
        std::sort(idx.begin(),idx.end(),[&](size_t aa,size_t bb){return e[aa]<e[bb];});

        hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);                // density matrix D = Sum_occ 2 c c^H (Hermitian)
        for (int k=0;k<Nocc;k++)
        {
            size_t c=idx[k];
            for (size_t i=0;i<n;i++)
                for (size_t j=i;j<n;j++)
                    D(i,j) += 2.0*dcmplx(U(i,c))*std::conj(dcmplx(U(j,c)));
        }
        Dprev = hmat_t<dcmplx>(0.6*Dprev + 0.4*D);               // linear mixing
        cd=std::make_unique<qchem::ChargeDensity::PeriodicIrrepCD<dcmplx>>(Dprev, &pw, irr);

        // GetTotalEnergy(new cd) computes the energy AND invalidates the dynamic terms' Irrep-keyed cache
        // (their GetEnergy calls newCD) so the NEXT GetMatrix rebuilds the Hartree/XC matrices fresh.
        E=ham.GetTotalEnergy(cd.get());
        double Etot=E.GetTotalEnergy();
        if (std::abs(Etot-Eprev)<1e-7) { converged=true; break; }
        Eprev=Etot;
    }
    ASSERT_TRUE(converged);
    std::cout << "[Si framework-Gamma] charge="<<cd->GetTotalCharge()<<" Etot="<<E.GetTotalEnergy()
              << "  (Ekin="<<E["Kinetic"]<<" Een="<<E["Een"]<<" Eee="<<E["Eee"]<<" Exc="<<E["Exc"]<<")" << std::endl;

    EXPECT_NEAR(cd->GetTotalCharge(), 8.0, 1e-6);                 // 8 valence electrons
    // This manual Hamiltonian has no ion-ion term; the band-structure (electronic) energy is the
    // prototype anchor (the dropped-G=0 alignment E_alphaZ is now separated out of it).
    EXPECT_NEAR(E.GetElectronicEnergy(), 1.468, 5e-3);           // matches the standalone prototype Si-Gamma
}


// The SAME Si-Gamma Kohn-Sham problem as FrameworkSiliconGammaMatchesPrototype, but now driven by the
// REAL framework cSCFIterator (no hand-rolled SCF loop): cSCFIterator -> cWaveFunction (the composite WF
// -> IrrepWF) -> SCFAcceleratorNull's <dcmplx> diagonalize -> TOrbitals<dcmplx> fill -> PeriodicIrrepCD<dcmplx>,
// with the cHamiltonianImp summing the PW terms.  This is the milestone that retires the "k-loop
// in the IBS": single-k plane-wave DFT IS now Hamiltonian = Sum terms + SCFIterator, like atoms/molecules.
TEST(PW_Si, Γ_Anchor)
{
    using namespace qchem::Hamiltonian;
    const double a=10.26;                       // Si conventional cubic lattice constant (a.u.)
    FCCUnitCell cell(a);                         // FCC primitive cell
    cell.AddAtom(14, {0,0,0});                   // Si diamond: itsZ = TRUE species Z = 14; Zion comes from the PP (ZionFn)
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D  lat(cell, ivec3_t(1,1,1));

    // Run the Si-Gamma SCF at plane-wave cutoff \a Ecut under a given SEED strategy, EVERYTHING ELSE
    // identical (config, Hamiltonian, complex-DIIS accelerator, convergence params) -- so the only difference
    // is the seed and the iteration counts are a fair head-to-head.  Returns iters/convergence/charge/energy.
    namespace L3=BasisSet::Lattice;
    struct Run { size_t iters; bool conv; double charge; qchem::EnergyBreakdown E; };
    auto run=[&](double Ecut, qchem::ChargeDensity::SeedStrategy seed)
    {
        using namespace qchem::Hamiltonian;
        // The basis is an abstract tBasisSet<dcmplx> owning its plane-wave Bloch block(s); the PP model lives
        // on the Hamiltonian term (the pseudo-wall), NOT the basis.
        std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, Ecut));
        Irrep      irr=bs->GetIrreps(Spin::None)[0];
        Crystal_EC ec(irr, 8);   // Si insulator: 8 valence electrons (2 atoms x Zion 4)
        // Plane-wave LDA Kohn-Sham Hamiltonian: kinetic + external(pseudo) + Hartree + Dirac X + VWN5
        // (heap; the SCFIterator takes ownership).  Atoms from the lattice; the PP model lives on the term.
        cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure(), bs.get(), "Si", "LDA", 4);   // one-call: looks up + owns the GTH PP
        std::unique_ptr<cHamiltonian> hamOwner(ham);   // R2.22: the iterator borrows; this scope owns
        // Complex DIIS (Pulay): extrapolate the Fock matrix from the [F,D] history; damps the marginal drift
        // within Si's degenerate Gamma_25' manifold so it converges to a tight |Delta rho|.
        using qchem::SCFAccelerators::DIISParams;
        auto* acc=new qchem::SCFAccelerators::SCFAcceleratorDIIS(DIISParams{8, 0.5, 1e-10, 1e-9});
        std::unique_ptr<qchem::SCFAccelerators::SCFAccelerator> accOwner(acc);   // R2.22: the iterator borrows; this scope owns
        qchem::SCFIterator::SolidSCFIterator scf(bs.get(), &ec, ham, acc, seed, lat.GetStructure().get());

        SCFParams par;
        par.NMaxIter      =80;
        par.MinΔρ         =1e-7;   // tight convergence (complex DIIS damps the degenerate-manifold drift)
        par.MinΔFD        =1e30;
        par.MinVirial     =1e30;   // pseudopotential calc: the textbook -V/K=2 virial does not hold
        par.MinFD         =1e30;
        par.StartingRelaxRo=0.4;
        par.MergeTol      =1e-4;
        par.Verbose       =false;
        scf.Iterate(par);

        auto cd=scf.GetWaveFunction()->GetChargeDensity();   // BUILT for us; the unique_ptr owns it
        double charge=cd->GetTotalCharge();
        return Run{scf.GetIterationCount(), scf.Converged(), charge, scf.GetEnergy()};
    };

    // The plane-wave SAD path (FourierSeedCD: a G-space form-factor sum of the atomic VALENCE densities from
    // atomic_valence_densities.json).  That file now holds the SMOOTH pseudo-valence Si density produced by
    // the Atom-PP (Ham_PP + KB nonlocal: scfrun --model PP --valence 4 --out ...), not the old all-electron
    // valence whose core peak injected spurious high-G content.  See doc/SCFSeedingPlan.md section 9.7.

    // Ecut=4 / Gamma (fast, near-jellium).  Both seeds; cross-check the SAD energy vs the standalone prototype.
    Run sad=run(4.0, qchem::ChargeDensity::SeedStrategy::SAD);
    Run uni=run(4.0, qchem::ChargeDensity::SeedStrategy::Uniform);
    std::cout << "[Si SCFIterator-Gamma] SAD iters="<<sad.iters<<"  Uniform iters="<<uni.iters
              << "  SAD Etot="<<sad.E.GetTotalEnergy()
              << "  (Ekin="<<sad.E["Kinetic"]<<" Een="<<sad.E["Een"]<<" Eee="<<sad.E["Eee"]<<" Exc="<<sad.E["Exc"]<<")" << std::endl;

    EXPECT_TRUE(sad.conv);
    EXPECT_NEAR(sad.charge,                  8.0,    1e-6);   // 8 valence electrons
    EXPECT_NEAR(sad.E.GetElectronicEnergy(), 1.468,  5e-3);   // band energy matches the standalone prototype
    // Physical total: electronic + ion-ion Ewald (Enn) + dropped-G=0 alignment (E_alphaZ) -- now NEGATIVE.
    // (Underconverged at Ecut=4 / Gamma-only; the converged Si total is ~-7.9 Ha/cell.)
    EXPECT_NEAR(sad.E.GetTotalEnergy(),     -7.2273, 5e-3) << "Enn="<<sad.E["Enn"]<<" E_alphaZ="<<sad.E["E_alphaZ"];
    EXPECT_TRUE(uni.conv);
    EXPECT_NEAR(sad.E.GetTotalEnergy(), uni.E.GetTotalEnergy(), 1e-3);   // seed cannot change the answer

    // THE SMOOTH-DENSITY WIN (regression guard).  The smooth pseudo-valence seed converges in ~12 iters here;
    // the OLD all-electron-valence seed needed ~15 (its core peak injected spurious high-G content).  Guard
    // the robust 3-iter improvement -- NOT a brittle tie with Uniform: at this coarse near-jellium Ecut=4 the
    // converged density is nearly flat so Uniform (the pure G=0 seed) is already near-optimal (11).  The
    // STRICT win over Uniform appears once the density has real structure -- e.g. Ecut=9/Gamma gives SAD 10
    // vs Uniform 11 (verified manually) -- but a 1-iter margin behind an ~18s high-Ecut SCF is too slow and
    // brittle to gate CI.  See doc/SCFSeedingPlan.md section 9.7.
    EXPECT_LE(sad.iters, 13u) << "smooth pseudo-valence SAD should converge in <=13 iters (was ~15 all-electron)";
}


// (The GTH database-reader unit test lives in src/BasisSet/Lattice/tests/GTH_UT.C -- it needs only
// the basis layer, not the SCF stack.)

// Silicon again, but BZ-sampled over a 2x2x2 Monkhorst-Pack mesh (8 k-points) instead of Gamma-only.
// This exercises the irrep(k) loop, the k-dependent KB projectors, and a properly BZ-averaged density.
// Validate the same robust internal properties (converges, exact 8-electron charge, band==direct
// energy) and report the BZ-sampled band gap.
TEST(PW_Si, k222_Prototype_Converges)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> A(0.0,h,h,  h,0.0,h,  h,h,0.0);
    UnitCell          cell(A);
    Lattice_3D        lat(cell, ivec3_t(2,2,2));
    ReciprocalLattice recip(lat.Reciprocal());
    const double      Omega=cell.GetCellVolume();
    const double      Ecut=4.0;

    Molecule si;
    si.Insert(new Atom(14, rvec3_t(0,0,0)));
    si.Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));
    GTH_PP                 siPP=GetGTH("Si","LDA",4);          // CP2K GTH-LDA q4, from the database
    const HGH_LocalPotential&     loc=siPP.local;
    const HGH_SeparablePotential& nl =siPP.nonlocal;

    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};

    // grid resolves the difference set: size from the richest (Gamma) basis.
    PlaneWave_IBS pwGamma(lat.Reciprocal(), ivec3_t(2,2,2), ivec3_t(0,0,0), Ecut);
    int m=MaxGComponent(pwGamma);
    ivec3_t Ng(4*m+1, 4*m+1, 4*m+1);

    SCFResult R=RunSCF_kpoints(recip, recip.GetCell(), Omega, si, loc, nl,
                               ivec3_t(2,2,2), Ecut, 8, Ng, vxcOf, epsOf, 0.4, 1e-8, 400);

    std::cout << "[Si 2x2x2] iters="<<R.iters<<" converged="<<R.converged
              << "\n  Etot(band)="<<R.Etot_band<<" Etot(direct)="<<R.Etot_direct
              << " gap="<<R.gap<<" Ha ("<<R.gap*27.2114<<" eV)" << std::endl;

    ASSERT_TRUE(R.converged);
    EXPECT_NEAR(Omega*std::real(RhoAt(R.rho,ivec3_t(0,0,0))), 8.0, 1e-6);   // 8 valence e-
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-5);                          // stationary fixed point
    EXPECT_GT (R.gap, 0.0);                                                 // Si is a semiconductor
}


// Stage 4 / multi-k: the SAME Si Kohn-Sham problem on a 2x2x2 Brillouin-zone mesh, through the REAL
// cSCFIterator.  Now the basis holds 8 Bloch blocks (one per k-point); the framework's per-irrep loop
// (MakeIrrepWFs, one IrrepWF per block) IS the BZ sum Sum_k w_k, with each block's density BZ-weighted
// (Symmetry::GetWeight = w_k = 1/8) so the total charge is 8 (not 8*8).  Reproduces the standalone
// prototype ScfSiliconBZSampled (Etot=0.934, gap>0).
TEST(PW_Si, k222_Anchor)
{
    using namespace qchem::Hamiltonian;
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});                   // itsZ = TRUE species Z = 14 (Zion=4 supplied by the PP's ZionFn)
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D  lat(cell, ivec3_t(2,2,2));       // 2x2x2 = 8 k-points

    namespace L3=BasisSet::Lattice;
    std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, 4.0));
    BasisSet::irrepv_t irreps=bs->GetIrreps(Spin::None);   // one Bloch irrep per k-block (8)

    const int Nelec=8;
    Crystal_EC ec(irreps, Nelec);                          // Nval per k-block; weights handle the BZ sum

    cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure(), bs.get(), "Si", "LDA", 4);   // one-call: looks up + owns the GTH PP
    std::unique_ptr<cHamiltonian> hamOwner(ham);   // R2.22: the iterator borrows; this scope owns
    using qchem::SCFAccelerators::DIISParams;
    auto* acc=new qchem::SCFAccelerators::SCFAcceleratorDIIS(DIISParams{8, 0.5, 1e-10, 1e-9});
    std::unique_ptr<qchem::SCFAccelerators::SCFAccelerator> accOwner(acc);   // R2.22: the iterator borrows; this scope owns

    // Uniform-density seed on the first block: D=(N/n0)I gives rho(r)=N/V (uniform), the total density
    // every block's first Hartree/XC needs (a single block suffices since rho is constant).  Built
    // centrally by MakeSeedDensity (Uniform); also the plane-wave Default.
    qchem::SCFIterator::SolidSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::Uniform);

    SCFParams par;
    par.NMaxIter      =80;
    par.MinΔρ         =1e-7;   // tight: the FFT XC/Hartree made the multi-k run cheap, so complex DIIS now
                               // converges the BZ-sampled density fully (Etot=0.934176 exactly) in seconds.
    par.MinΔFD        =1e30;
    par.MinVirial     =1e30;
    par.MinFD         =1e30;
    par.StartingRelaxRo=0.4;
    par.MergeTol      =1e-4;
    par.Verbose       =false;
    scf.Iterate(par);

    const qchem::WaveFunction::cWaveFunction* wf=scf.GetWaveFunction();
    auto cd=wf->GetChargeDensity();
    double charge=cd->GetTotalCharge();
    qchem::EnergyBreakdown E=scf.GetEnergy();

    // Gap = lowest unoccupied - highest occupied, merged across all k-blocks (occupations are physical:
    // 2 for filled bands, 0 above, since only the density is BZ-weighted, not the occupation).
    double homo=-1e30, lumo=1e30;
    for (const auto& [e,lev] : wf->GetEnergyLevels())
        if (lev.occ>1e-6) homo=std::max(homo,e); else lumo=std::min(lumo,e);
    double gap=lumo-homo;

    std::cout << "[Si SCFIterator-2x2x2] iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Etot="<<E.GetTotalEnergy()<<" gap="<<gap<<" Ha ("<<gap*27.2114<<" eV)"
              << "  (Ekin="<<E["Kinetic"]<<" Een="<<E["Een"]<<" Eee="<<E["Eee"]<<" Exc="<<E["Exc"]<<")" << std::endl;

    EXPECT_TRUE(scf.Converged());
    EXPECT_NEAR(charge,                  8.0,    1e-6);    // 8 valence electrons (BZ-weighted sum)
    EXPECT_NEAR(E.GetElectronicEnergy(), 0.934,  5e-3);    // band energy matches prototype ScfSiliconBZSampled
    // Physical total = electronic + ion-ion Ewald (same per-cell Madelung as Gamma) + G=0 alignment.
    EXPECT_NEAR(E.GetTotalEnergy(),     -7.7613, 5e-3) << "Enn="<<E["Enn"]<<" E_alphaZ="<<E["E_alphaZ"];
    EXPECT_GT(gap, 0.0);                              // Si is a semiconductor
}

} //namespace
