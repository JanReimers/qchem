// File: Calculation/Imp/MoleculeCalculation.C  Implementation of the qchem::MoleculeCalculation facade.
module;
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <nlohmann/json.hpp>
module qchem.MoleculeCalculation;

import qchem.BasisSet.Gaussian.Point.Factory;             // BasisSet::Gaussian::Factory
import qchem.BasisSet.Gaussian.Point.BasisFiles;          // BasisFileChecksum: a stated basis file's identity (report, fingerprint)
import qchem.BasisSet.Conditioning;                  // ReportConditioning
import qchem.BasisSet.IntegralsCacheReport;          // EmitIntegralsCacheReport (the run's `cache` section)
import qchem.BasisSet.Gaussian.PG_Cart.SymmetryAdapt;   // PG::SymmetryAdapt (the SALC builder)
import qchem.ElectronConfiguration.Molecule;        // Molecule_EC (global aufbau)
import qchem.SCFAccelerator.Factory;                // SCFAccelerators::Type, Factory
import qchem.WaveFunction;                           // WaveFunction (GetChargeDensity/GetOrbitals/GetQNs)
import qchem.Orbitals;                               // Orbital, Orbitals
import qchem.Reporting;                              // report:: run sink + the scf-section writer
import qchem.SCFRecipeReport;                        // ReportSCFRecipe: the `scf` header every facade files (V-FACADE-PARITY)
import qchem.DegenerateFillReport;                   // the `choose` row for a kT = 0 fill that shared a degenerate level
import qchem.BasisSet.Orbital_1E_IBS;                // Real_OIBS -- the per-irrep blocks a carry is re-homed onto
import qchem.CompositeCD;                            // tComposite_CD (the carried density's shape)
import qchem.ChargeDensity.Factory;                  // IrrepCD_Factory
import qchem.Symmetry.Molecule.PointGroup;           // DetectPointGroup (a session step must keep the SALC group)
import qchem.Hamiltonian.Factory;
import qchem.Response;          // MakeReference, MakeOrbitalFrame, MakeDipoleProbe
import qchem.Response.Solver;   // LinearResponse
import qchem.BasisSet.AoShellSource;   // MeshAudit: the largest element named by atom, shell and l
import qchem.Mesh.BeckeClimb;          // ClimbBeckeLadder: the mesh pilot's climb at a fixed D
import qchem.Policy;                   // policy::Report: the pilot's decision files itself

namespace qchem
{

namespace PG = ::qchem::BasisSet::Gaussian::PG_Cart;
using SCFIter = qchem::SCFIterator::MolecularSCFIterator;   // item 2: the atom/molecule display subclass

//! The spin density of an unpolarized run: m(r) = 0 everywhere (GetSpinDensity).
struct ZeroField : public virtual ScalarFunction<double>
{
    double         operator()(const rvec3_t&) const override {return 0.0;}
    vec3_t<double> Gradient  (const rvec3_t&) const override {return vec3_t<double>(0.0, 0.0, 0.0);}
};

// Build the molecular orbital basis, optionally SALC point-group-blocked.  When symmetry is off the
// caller owns the raw basis directly; when on, the raw basis is handed to PG::SymmetryAdapt as a
// shared_ptr that the returned SymmetryAdaptedBasisSet keep-alives -- so the facade owns (and deletes)
// only the adapted basis, and the raw lifetime rides along.  PG::SymmetryAdapt throws if the basis has
// no Cartesian PGData IBS (the guard); the facade only ever builds that basis today, so it is safe.
// \a ao receives the SALC run's raw basis (shared with the adapted one), whose shells identify a saved state; null without symmetry.
static BasisSet::Real_BS* BuildBasis(const CalcOptions& opts, const std::shared_ptr<const Structure>& st,
                                     std::shared_ptr<const BasisSet::Real_BS>& ao)
{
    ao.reset();
    const bool spherical = (opts.angular == Angular::Spherical);
    // SALC now supports the Cartesian PGData basis AND the in-house MnD-spherical basis (SphData;
    // doc/OldPlans/SphericalSALCPlan.md S3a/S4).  Only libcint-spherical is unwired (S3b) -- and it presents as a
    // PGData with spherical components, which the Cartesian extractor would silently misread, so reject that
    // one combination clearly here rather than let SymmetryAdapt take the wrong path.
    if (opts.symmetry && spherical && opts.engine == Engine::LibCint)
        throw std::runtime_error("qchem::MoleculeCalculation: {.symmetry=true} with angular=Spherical requires "
                                 "engine=MnD (libcint-spherical SALC is not yet supported)");

    const char* engine  = (opts.engine  == Engine::LibCint) ? "libcint" : "mnd";
    const char* angular = spherical ? "spherical" : "cartesian";
    BasisSet::Real_BS* raw = opts.basisFile.empty()
        ? BasisSet::Gaussian::Factory(nlohmann::json{{"basis", opts.basis}, {"engine", engine}, {"angular", angular}}, st.get())
        : BasisSet::Gaussian::FactoryFromFile(opts.basisFile, st.get(),                               // a stated file (CalcOptions::basisFile)
              opts.engine==Engine::LibCint ? BasisSet::Gaussian::Engine::LibCint : BasisSet::Gaussian::Engine::MnD,
              spherical ? BasisSet::Gaussian::Angular::Spherical : BasisSet::Gaussian::Angular::Cartesian);
    if (!opts.symmetry) return raw;
    std::shared_ptr<const BasisSet::Real_BS> rawShared(raw);
    ao=rawShared;
    return PG::SymmetryAdapt(rawShared, *st, opts.symmetryTol);
}

// Build the molecular electron configuration from the requested multiplicity 2S+1 (validated by MoleculeSetup::Resolve).
// multiplicity<=0 is the minimal-spin default (closed-shell singlet for even Ne, doublet for odd, via Molecule_EC(Ne)).
static Molecule_EC* MakeMoleculeEC(int Ne, int multiplicity)
{
    if (multiplicity <= 0) return new Molecule_EC(Ne);    // auto / minimal spin (historical behaviour)
    const int twoS = multiplicity - 1;                    // = nUp - nDown
    return new Molecule_EC((Ne + twoS) / 2, (Ne - twoS) / 2);   // (nUp, nDown)
}

MoleculeCalculation::MoleculeCalculation(std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf, Observer onIteration)
    : MoleculeCalculation(setup, *setup->reference, nullptr, scf, std::move(onIteration))
{}

MoleculeCalculation::MoleculeCalculation(std::shared_ptr<const MoleculeSetup> setup, const Structure& at, std::unique_ptr<const Carry> carry,
                         const SCFParams& scf, Observer onIteration, bool geometryPoint)
    : MoleculeCalculation(BuildOnly{}, std::move(setup), at, std::move(carry), std::move(onIteration), geometryPoint)
{
    (void)Converge(scf);                              // so Energy()/Density() are ready on return; Result() is its outcome
}

std::unique_ptr<MoleculeCalculation> MoleculeCalculation::Prepare(std::shared_ptr<const MoleculeSetup> setup, Observer onIteration)
{
    const Structure& ref=*setup->reference;
    std::unique_ptr<MoleculeCalculation> c(new MoleculeCalculation(BuildOnly{}, std::move(setup), ref, nullptr, std::move(onIteration)));
    c->BuildGraph();                                  // into the CALLER'S run (user ruling D, 2026-11-02: drivers own the runs)
    c->itsPrepared=true;
    return c;
}

MoleculeCalculation::MoleculeCalculation(BuildOnly, std::shared_ptr<const MoleculeSetup> setup, const Structure& at,
                         std::unique_ptr<const Carry> carry, Observer onIteration, bool geometryPoint)
    : itsGeometryPoint(geometryPoint)
    , itsSetup(std::move(setup))
    // A pseudopotential run works on the VALENCE ions (Z-Zion charge) so the electron count is the valence count -- the
    // structure the basis and the EC both read.
    , itsStructure(itsSetup->PointStructure(at))
    , itsCarry(std::move(carry))
    , itsObserver(std::move(onIteration))
{
    itsEC    = MakeMoleculeEC(itsSetup->Ne, itsSetup->opts.multiplicity);
    // The basis is NOT built here: BuildGraph builds it INSIDE an open run report's "basis" scope, so a basis with something
    // to announce does so as it is constructed (V1.14).
}

MoleculeCalculation::~MoleculeCalculation()
{
    delete itsScf;     // FIRST: the iterator references itsHam/itsAccel, which are freed after this body (R2.22)
    delete itsBasis;
    delete itsEC;
}

// THE GRAPH, built inside the caller's open report run (if any: with none, every report write is a no-op).
void MoleculeCalculation::BuildGraph()
{
    namespace rpt = qchem::report;
    const CalcOptions& opts=itsSetup->opts;
    const nlohmann::json jsacc = ToAcceleratorJson(itsSetup->acc, itsSetup->eMax);   // every knob of the chain, forwarded
    // The facade sets the context, the accelerator files what it was built with (THE RUN REPORT): one `accelerators` row per
    // accelerator built in this run (a solid rebuilds one per stage), its Ladder rungs one level further in.
    nlohmann::json jsrun = jsacc;
    if (itsStartOnGDM) { jsrun["startOnGDM"]=true; itsStartOnGDM=false; }   // the restored stage only (PrepareRestart)
    qchem::SCFAccelerators::SCFAccelerator* accel=nullptr;
    { rpt::Row built("accelerators"); accel = qchem::SCFAccelerators::Factory(itsSetup->accelerator, jsrun); }
    const auto seed = itsSetup->seed;

    // The `basis` section is ASSEMBLED across layers, and since V1.14 it spans the basis's whole LIFE in this
    // run: the FIRST Converge BUILDS the basis (a basis with something to announce does so as it is constructed --
    // nobody asks it to), the orchestrator stamps the scalars and asks for the conditioning
    // report (basis.perIrrep, removed, conditioning).  The HAMILTONIAN is built OUTSIDE it, in its own context "hamiltonian" (user rulings 2026-11-02:
    // the caller sets the context, the callee reports one level in; 2026-11-04: the report's nesting is the run sequence)
    // -- each mesh files there as it is built, just before the term that needs it, and each term its fit basis, so a
    // molecule and a solid file the same paths.  So the basis context is entered twice: SILENTLY for the build (a Scope renders only the
    // scalars it added, and the build adds rows), then as the rendering Section for the rest -- which renders the
    // whole basis once, complete, on close.  A later Converge reuses the basis.
    {
        rpt::Scope basis("basis");
        if (!itsBasis) itsBasis = BuildBasis(opts, itsStructure, itsAoBasis);   // raw, or SALC-blocked if .symmetry
    }
    // The unified resolver turns the Model token into the concrete Hamiltonian; the DFT extras (mesh, orbital basis,
    // xalpha) are ignored for HF/1-e/Dirac.  A pseudopotential run takes the PP front door instead (LSDA valence
    // Hamiltonian: V_loc + KB projectors + Zion ion-ion in place of Ven).
    decltype(itsHam.get()) ham=nullptr;
    {
        rpt::Scope hamiltonian("hamiltonian");
        ham = MakeHamiltonianOn(itsSetup->mesh);
    }
    {
        rpt::Section basis("basis");
        rpt::Set("name",       opts.basis);
        if (!opts.basisFile.empty())
        {
            rpt::Set("file",     opts.basisFile);
            rpt::Set("checksum", BasisSet::Gaussian::BasisFileChecksum(opts.basisFile));
        }
        rpt::Set("engine",     opts.engine  == Engine::LibCint  ? "libcint"   : "mnd");
        rpt::Set("angular",    opts.angular == Angular::Spherical ? "spherical" : "cartesian");
        rpt::Set("nFunctions", (long)itsBasis->GetNumFunctions());
        qchem::BasisSet::ReportConditioning(*itsBasis);
        // R2.22: the iterator no longer deletes what it is handed, so this facade adopts the pair.  Order is
        // deliberate -- the PREVIOUS iterator dies first, then the previous Hamiltonian/accelerator it pointed
        // at are freed by these resets, which is exactly when `delete itsScf` used to free them.
        delete itsScf; itsScf = nullptr;
        itsHam.reset(ham);
        itsAccel.reset(accel);
        if (itsCarry && itsCarry->moved) itsCarry=std::make_unique<const Carry>(CarryOntoThisBasis(*itsCarry));   // a session step
        if (!itsCarry)
            itsScf = new SCFIter(itsBasis, itsEC, itsHam.get(), itsAccel.get(), seed, itsStructure.get());
        else
        {
            // THE CARRY (a MoleculeSession step, or a restored state): the D, block for block, onto THIS point's basis -- a
            // WorkingDensity of the moved Hamiltonian, or an EXACT RESUME of a state saved by this very run, consumed here once.
            itsScf = new SCFIter(itsBasis, itsEC, itsHam.get(), itsAccel.get(), MakeCarriedDensity(*itsCarry).release(), itsStructure.get(),
                                 qchem::Auto, 0.0, itsCarryExact ? qchem::SCFIterator::SeedRole::ExactResume
                                                                 : qchem::SCFIterator::SeedRole::WorkingDensity);
            itsCarry.reset();
            itsCarryExact=false;
        }
    }
}

Outcome<MoleculeCalculation::Converged, SCFFailure> MoleculeCalculation::Converge(const SCFParams& asked)
{
    // The SCF writes into the CALLER'S report run: a facade never opens one (user ruling D, 2026-11-02 -- drivers own the
    // runs, the caller sets the context one level up).  With no run open every write is a no-op.
    namespace rpt = qchem::report;
    // THE DEGENERATE-LEVEL RULE (D-DEGEN-FILL, user 2026-11-08): a stated rule passes; Auto is THIS facade's policy, by model.
    //   * DFT (LDA, Xalpha, any PP run): Shared, the kT -> 0 Fermi limit (D invariant to the rotation inside the level);
    //   * HF and the one-electron models: OrbitalByOrbital (user ruling 2026-11-08, "if B gives the lower symmetry, lower energy, then
    //     B"): HF charges a fractional open shell its same-spin self-interaction, so the integer symmetry-broken state lies lower
    //     (Boron UHF/dzvp: -24.525932 integer vs -24.40769 shared).  Which of the equivalent orbitals fills is arbitrary among them.
    using qchem::Hamiltonian::Model;
    const Model model=itsSetup->opts.model;
    const bool isDFT = itsSetup->pp || model==Model::LDA || model==Model::Xalpha;
    const SCFParams params=ResolveDegenerateLevelFill(asked,
        isDFT ? SCFParams::DegenerateLevelFill::Shared : SCFParams::DegenerateLevelFill::OrbitalByOrbital,
        isDFT ? "a DFT model: a partly filled degenerate level takes the kT -> 0 Fermi limit, so its density does not depend on the "
                "eigensolver's rotation inside the level"
              : "HF: integer fill, the lower-energy symmetry-broken state; which degenerate orbital fills is arbitrary among equivalent ones "
                "(a fractionally shared HF open shell pays its same-spin self-interaction: Boron UHF/dzvp 118 mHa higher).  State "
                "degenerateLevelFill: Shared for the symmetric state", true,
        // A GDM run (user 2026-11-08: guide the user; GDM works only on integer occupations): a molecule mixes in D space.
        DegenerateFillContext{.directMin = itsSetup->accelerator==SCFAccelerators::Type::GDM,
                              .fallback  = DescribeDeclinedFallback(ResolveMixing(asked, /*gSpace*/false))});
    {   // the "what recipe am I running" header; eps is consumed by the XC mesh wherever there is one (pin 31)
        rpt::Section scf("scf"); ReportSCFRecipe(params, itsSetup->acc.type);
        ReportAccuracy(itsSetup->opts.accuracy, !HasXCMesh()           ? std::string(kNoAccuracyConsumerYet)
                                              : itsSetup->meshResolved ? std::string("the Becke XC mesh (sized by an earlier pilot)")
                                              : itsSetup->meshStated   ? std::string("the Becke XC mesh (stated: audited, never raised)")
                                              :                          std::string("the Becke XC mesh (climbed at the converged density)"));
    }
    if (itsPrepared) itsPrepared=false;                       // Prepare built the graph (in its caller's run): reuse it, once
    else             BuildGraph();
    itsLastkT=params.SmearingkT;
    AttachObserver();
    (void)itsScf->Iterate(params);   // the outcome is the iterator's Converged(), minted below
    RebuildSampling();
    if (!itsPiloted && HasXCMesh() && !itsSetup->meshResolved && itsScf->Converged()) PilotMesh(params);

    // Snapshot the integrals cache into the run's `cache` section (RunReportPlan step 4), inside the bracket.
    // The cache is a process-wide singleton NEVER cleared between runs, so this is cumulative-to-this-point --
    // exact for a one-run process, the running total in a multi-run one.
    qchem::BasisSet::EmitIntegralsCacheReport();

    ReportResult("molecule", itsScf->Converged(), itsScf->GetIterationCount(), itsDensity->GetTotalCharge(), itsScf->GetEnergy(),
                 itsSetup->opts.accuracy.Target(Consumes(), itsStructure->GetNumAtoms()));
    return Outcome_();
}

Outcome<MoleculeCalculation::Converged, SCFFailure> MoleculeCalculation::Result() const {return Outcome_();}

Outcome<MoleculeCalculation::Converged, SCFFailure> MoleculeCalculation::Outcome_() const
{
    using O=Outcome<Converged, SCFFailure>;
    if (itsScf->Converged()) return O::Ok(Converged(this));
    SCFFailure f;   // the molecular facade's detectors: the iterator's own refusal, else a run out of iterations
    const bool refused=!itsScf->UnphysicalStop().empty();
    f.why        = refused ? SCFFailure::Why::UnphysicalEnergy : SCFFailure::Why::NotConverged;
    f.iterations = itsScf->GetIterationCount();
    f.lastEnergy = itsScf->GetEnergy().GetTotalEnergy();
    f.details    = refused ? itsScf->UnphysicalStop() : "the SCF did not converge in "+std::to_string(f.iterations)+" iterations";
    if (!refused && itsLastRow) f.details += "; the last iterate: "+itsLastRow->GatesText();
    return O::Fail(std::move(f));
}

double                 MoleculeCalculation::Converged::Energy()         const {return itsCalc->itsScf->GetEnergy().GetTotalEnergy();}
qchem::EnergyBreakdown MoleculeCalculation::Converged::EnergyTerms()    const {return itsCalc->itsScf->GetEnergy();}
size_t                 MoleculeCalculation::Converged::IterationCount() const {return itsCalc->itsScf->GetIterationCount();}
const MoleculeCalculation::sf_t& MoleculeCalculation::Converged::Density()      const {return itsCalc->Density();}

// The linear response of the converged state: Reference (R0, orbital basis) + frame (AO <-> MO) + the
// Hamiltonian's kernel + the dipole probe, solved together by GMRES (qcResponse).  The facade is the one place
// that owns all four ingredients, which is why the entry point is here.
Outcome<rmat_t,Response::ResponseFailure> MoleculeCalculation::StaticPolarizability(const qcMesh::MeshParams& dipoleMesh,
                                                                            const KrylovParams& krylov) const
{
    const auto* wf=itsScf ? itsScf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("MoleculeCalculation::StaticPolarizability: no converged wave function");
    // Integer occupations (the molecular aufbau fills across irreps, per spin channel).  The eigenvalue noise is
    // not measured by this facade yet: NaN says so, and E1 then gates on the gap's sign only.
    const Response::Reference ref=Response::MakeReference(*wf, OccupationConfig{}, {.acrossK=true, .acrossSpin=false},
                                                          std::numeric_limits<double>::quiet_NaN());
    const Response::OrbitalFrame<double> frame=Response::MakeOrbitalFrame(ref, *wf);
    const auto kernel=itsHam->MakeResponseKernel(itsBasis, itsDensity.get());
    const qcMesh::Mesh probeMesh=[&]{ qchem::report::Scope here("dipoleProbe");   // the mesh builder reports one level in
                                      return itsStructure->CreateIntegrationMesh(dipoleMesh); }();
    const Response::OperatorProbe probe=Response::MakeDipoleProbe(ref, frame, *wf, probeMesh);
    auto r=Response::LinearResponse(ref, frame, *kernel, probe, std::make_shared<Symmetry::Invariant>(), krylov);
    using O=Outcome<rmat_t,Response::ResponseFailure>;
    if (!r)
    {
        std::cout << "[response] FAILED: " << r.Error().detail << std::endl;
        return O::Fail(r.Error());
    }
    r->Write(std::cout);
    // H = H0 + F.r for an electron (charge -1), so mu = -<r> and alpha = d mu / dF = -chi.
    rmat_t alpha(3,3);
    for (size_t i=0;i<3;i++) for (size_t j=0;j<3;j++) alpha(i,j)=-r->chi(i,j).real();
    return O::Ok(std::move(alpha));
}

void MoleculeCalculation::RebuildSampling()
{
    const auto* wf = itsScf->GetWaveFunction();
    itsDensity = wf->GetChargeDensity();               // it BUILDS the density and hands ownership over

    itsOccupied.clear();
    for (const Irrep& irr : wf->GetQNs())
        for (const auto* o : wf->GetOrbitals(irr)->Iterate())   // Iterate() yields const Orbital*
            if (o->IsOccupied())
            {
                // An Orbital IS a ScalarFunction<double> (cast between abstract faces, per CLAUDE.md).
                const auto* phi = dynamic_cast<const sf_t*>(o);
                assert(phi);
                itsOccupied.push_back({o->GetEigenEnergy(), phi});
            }
    std::sort(itsOccupied.begin(), itsOccupied.end(),
              [](const occ_t& a, const occ_t& b){ return a.first < b.first; });   // ascending energy

    itsLevels.clear();
    for (const Irrep& irr : wf->GetQNs())
        for (const auto* o : wf->GetOrbitals(irr)->Iterate())
        {
            const auto* phi = dynamic_cast<const sf_t*>(o);   // abstract -> abstract, as above
            assert(phi);
            itsLevels.push_back({o->GetEigenEnergy(), o->GetOccupation(), irr.ms, o->GetLabel(), phi});
        }
    std::stable_sort(itsLevels.begin(), itsLevels.end(), [](const Level& a, const Level& b){ return a.e < b.e; });

    if (wf->GetSpinGroup()==SpinGroup::Polarized)
        itsSpinDensity = wf->GetSpinDensity();
    else
        itsSpinDensity = std::make_unique<ZeroField>();   // m == 0 under imposed SU(2): a field, not an absence
}

void MoleculeCalculation::OnIteration(Observer obs) {itsObserver = std::move(obs);}

std::vector<hmat_t<double>> MoleculeCalculation::BlockDensities(bool energyWeighted) const
{
    const auto* wf = itsScf ? itsScf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("qchem::MoleculeCalculation: no wave function yet (a block density needs a Converge)");
    std::vector<hmat_t<double>> out;
    for (Spin s : SpinIrreps(itsSetup->opts.spin))
        for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            const size_t n=ibs->GetNumFunctions();
            const auto* os=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(wf->GetOrbitals(ibs->GetIrrep(s)));
            if (!os) throw std::logic_error("qchem::MoleculeCalculation: a block's orbitals are not real");
            hmat_t<double> D(n);
            for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) D(a,b)=0.0;
            for (const auto* o : os->Iterate<qchem::Orbitals::TOrbital<double>>())
            {
                const double f=o->GetOccupation()*(energyWeighted ? o->GetEigenEnergy() : 1.0);
                if (o->GetOccupation()==0.0) continue;
                const auto& v=o->GetCoeff();
                for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) D(a,b)+=f*v[a]*v[b];
            }
            out.push_back(std::move(D));
        }
    return out;
}

MoleculeCalculation::Carry MoleculeCalculation::MakeCarry() const {return Carry{BlockDensities(false)};}
MoleculeCalculation::Carry MoleculeCalculation::MakeMovedCarry(double displacement) const
{
    return Carry{.D={}, .moved=SnapshotState(StateBlocks(true), StateKind()), .displacement=displacement};
}

MoleculeCalculation::Carry MoleculeCalculation::CarryOntoThisBasis(const Carry& carry)
{
    assert(carry.moved);
    SavedState st=*carry.moved;
    const std::vector<StateBlock> now=StateBlocks(false);
    CarryToMovedAtoms(st, now);   // the AO coefficients follow their atoms (the shared carry; its declaration has the measurement)
    auto D=CheckedBlockDensities(st, now);
    // Both points come from ONE setup (same recipe, same blocks), so a refusal here is a broken invariant, not a caller's choice.
    if (!D) throw std::logic_error("qchem::MoleculeCalculation: a session step's carry does not fit the next point: "+D.Error().details);
    Carry out;
    for (const mat_t<dcmplx>& M : *D)
    {
        const size_t n=M.rows();
        hmat_t<double> h(n);
        for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) h(a,b)=M(a,b).real();
        out.D.push_back(std::move(h));
    }
    std::ostringstream moved;
    moved<<"structure.positions: the atoms moved by up to "<<std::setprecision(6)<<carry.displacement<<" bohr from the previous point";
    itsRestartedFrom=RestartSummary{.path="(previous geometry point)", .exact=false, .warm={moved.str()},
                                    .maxDisplacement=carry.displacement, .electrons=Electrons(st, now)};
    report::FormatGuard keep(std::cout);
    std::cout<<"[molecule session] the previous point's density is carried as AO coefficients (atoms moved up to "
             <<std::setprecision(6)<<carry.displacement<<" bohr): "<<std::setprecision(12)<<itsRestartedFrom->electrons
             <<" electrons on this basis"<<std::endl;
    return out;
}
std::vector<hmat_t<double>> MoleculeCalculation::MakeEnergyWeightedDensity() const {return BlockDensities(true);}

std::unique_ptr<qchem::ChargeDensity::rDM_CD> MoleculeCalculation::MakeCarriedDensity(const Carry& carry) const
{
    using namespace qchem::ChargeDensity;
    auto cd=std::make_unique<tComposite_CD<double>>(itsBasis->GetImposition() ? itsBasis->GetImposition()->DensityGroup() : nullptr);
    size_t b=0;
    for (Spin s : SpinIrreps(itsSetup->opts.spin))
        for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            if (b>=carry.D.size() || carry.D[b].rows()!=ibs->GetNumFunctions())
                throw std::logic_error("qchem::MoleculeCalculation: the carried density's blocks do not match this point's basis");
            const Irrep irr=ibs->GetIrrep(s);
            cd->Insert(std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(carry.D[b++], ibs, irr)), irr);
        }
    if (b!=carry.D.size())
        throw std::logic_error("qchem::MoleculeCalculation: the carried density has more blocks than this point's basis");
    return cd;
}

GridAuditResult MoleculeCalculation::MeshAudit(const qcMesh::MeshParams& reference) const
{
    namespace rpt = qchem::report;
    const auto* wf = itsScf ? itsScf->GetWaveFunction() : nullptr;
    if (!wf || !itsDensity) throw std::logic_error("qchem::MoleculeCalculation::MeshAudit: no converged density yet (the audit runs at a fixed D)");
    const CalcOptions& opts=itsSetup->opts;
    const qchem::ChargeDensity::rDM_CD* d=itsDensity.get();
    const std::vector<hmat_t<double>> D=BlockDensities(false);
    GridAuditResult R;
    std::vector<std::vector<Symmetry::Molecule::AoShell>> shellsOf;   // per basis block: its AO shells, empty for a SALC block
    R.reference = "Becke mesh ("+std::to_string(reference.nRadial)+", "+std::to_string(reference.angularDegree)+")";
    {
        // FENCED: the reference Hamiltonian files its meshes and fit bases as it builds and first evaluates; a nested run keeps them out
        // of this run's `hamiltonian` section.  RAII: an exception mid-audit must not leave the report nested.
        struct Fence { explicit Fence(const std::string& n) {rpt::Begin(n);} ~Fence() {rpt::End();} } fence("meshAudit reference ("+R.reference+")");
        std::unique_ptr<qchem::Hamiltonian::rHamiltonian> hamRef(MakeHamiltonianOn(reference));
        itsHam->RefreshForDensity(itsBasis, d);
        hamRef->RefreshForDensity(itsBasis, d);
        const std::vector<qchem::Hamiltonian::TermEnergy> eRun=itsHam->MakeTermEnergies(d), eRef=hamRef->MakeTermEnergies(d);
        if (eRun.size()!=eRef.size()) throw std::logic_error("MeshAudit: the reference Hamiltonian has a different term list");
        for (size_t t=0;t<eRun.size();t++)
        {
            // Matched by the ENERGY ENTRIES, not the label: a mesh term's description states its own point count, which is what differs.
            if (eRun[t].entries!=eRef[t].entries)
                throw std::logic_error("MeshAudit: term "+std::to_string(t)+" reports '"+eRun[t].entries+"' here, '"+eRef[t].entries+"' on the reference");
            R.terms.push_back({eRun[t].label, eRun[t].entries, 0.0, eRun[t].E-eRef[t].E, eRun[t].E});
            R.ΔE+=eRun[t].E-eRef[t].E;
        }
        size_t channel=0;   // the carry's block order: spin irreps outer, basis blocks inner (BlockDensities')
        for (Spin s : SpinIrreps(opts.spin))
        {
            size_t block=0;
            for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
            {
                const std::unique_ptr<LASolver<double>> frame(LASolver<double>::Factory(qchem::Cholesky));
                frame->SetBasisOverlap(ibs->Overlap());
                const std::string spin = s==Spin::Up ? "up" : s==Spin::Down ? "down" : "none";
                FoldΔH(itsHam->MakeTermMatrices(ibs,s,d), hamRef->MakeTermMatrices(ibs,s,d), D.at(channel), 1.0, *frame, spin, block, R.terms);
                if (shellsOf.size()<=block)   // a block's shells, once (the spin channels share the blocks); a SALC block is not an AO block
                {
                    const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(ibs);
                    shellsOf.push_back(src && !opts.symmetry ? src->GetAoShells() : std::vector<Symmetry::Molecule::AoShell>{});
                }
                ++block; ++channel;
            }
        }
    }
    for (auto& t : R.terms) if (t.maxAbsΔH>0.0 && !shellsOf.at(t.maxElement.block).empty())   // the largest element of each term, named
    {
        t.maxElement.functionI=FunctionName(shellsOf[t.maxElement.block], *itsStructure, t.maxElement.i);
        t.maxElement.functionJ=FunctionName(shellsOf[t.maxElement.block], *itsStructure, t.maxElement.j);
    }
    for (auto& t : R.terms)
    {
        R.maxAbsΔH=std::max(R.maxAbsΔH, t.maxAbsΔH);
        if (t.entries=="Exc" && t.maxAbsΔH>0.0) t.etaOn="measured";   // the XC mesh has no rigorous bound
    }
    const ClientConsumes c=Consumes();
    ReportGridAudit(R, "[mesh audit]", "Becke mesh ("+std::to_string(itsSetup->mesh.nRadial)+", "+std::to_string(itsSetup->mesh.angularDegree)+")",
                    c==ClientConsumes::Energy ? 'E' : 'H', opts.accuracy.Target(c, itsStructure->GetNumAtoms()));
    return R;
}

qchem::Hamiltonian::rHamiltonian* MoleculeCalculation::MakeHamiltonianOn(const qcMesh::MeshParams& mesh) const
{
    const CalcOptions& opts=itsSetup->opts;
    return itsSetup->pp
        ? qchem::Hamiltonian::Factory(opts.spin, itsStructure, *itsSetup->pp, mesh, itsBasis, opts.vxcFit, opts.coulomb)
        : qchem::Hamiltonian::Factory(opts.model, opts.spin, itsStructure, mesh, itsBasis, opts.xalpha, opts.vxcFit, opts.coulomb);
}

bool MoleculeCalculation::HasXCMesh() const {return qchem::Hamiltonian::IsDFT(itsSetup->opts.model) || itsSetup->pp.has_value();}

ClientConsumes MoleculeCalculation::Consumes() const {return itsGeometryPoint ? ClientConsumes::MatrixElements : ClientConsumes::Energy;}

void MoleculeCalculation::AttachObserver()
{
    itsLastRow.reset();
    itsScf->SetObserver([this](const qchem::SCFIterator::IterationRow& row)   // keep the last row: a failure states its gates
    {
        itsLastRow = row;
        if (itsObserver) itsObserver(row);
    });
}

namespace
{
std::string MeshName(const qcMesh::MeshParams& m)
{
    return "("+std::to_string(m.nRadial)+", "+std::to_string(m.angularDegree)+(m.angular==qcMesh::AngularKind::GaussLegendre ? " GL" : "")+")";
}
//! Does the Hartree term's error FALL as the mesh rises, so a climb can size it?  The δ-ρ direct sum's does not: a Becke mesh is a poor rule
//! for the singular kernel near each point, and its J error sits at about 1e-2 on every rung tried (TeF₆: no better at 80 radial points than
//! at 50; Fitting/Internal/DirectSumCoulomb.C).  Chasing it raised water to the ladder's top at O(N²) per reference and per re-convergence
//! iteration (lane B, 2026-11-08: M_DFT.WaterLDA_DeltaCoulomb_Long ran ~10 h where its entry rung takes 196 s).  Every other route's
//! mesh-borne Hartree (δ-ρ by radial Poisson, the multipole channels) converges with the mesh and is climbed.
bool HartreeConvergesWithMesh(const Hamiltonian::CoulombRoute& c)
{
    return !(c.fit==Hamiltonian::CoulombFit::Delta && c.solver==Hamiltonian::CoulombSolver::DirectSum);
}
}

// THE MESH PILOT (doc/GeometryOptimizationPlan.md §14.2; pin 31).  The run has converged at its entry rung.  Climb the Becke ladder AT
// THIS D (qchem.Mesh.BeckeClimb, the solid's rule): each step compares the rung against the rung above on each axis, a molecule's basis
// does not depend on its mesh, so both sides are the SAME basis and the SAME D, and nothing is restored.  The distance is the run's own
// measure (Consumes), the worst over the mesh-borne terms (XC, and a δ-ρ Hartree rung riding the mesh; every other term reads 0), less a
// Hartree whose error does not fall with the mesh (HartreeConvergesWithMesh: the direct sum, WARNed, never chased).  If the
// climb raised: rebuild the point on the settled rung, re-converge ONCE from the carried D, and audit again (one raise at most) -- up to
// kPilotRounds.  A STATED mesh is audited once and never raised.  Then FREEZE: SizedSetup() carries the settled rung, resolved.
void MoleculeCalculation::PilotMesh(const SCFParams& params)
{
    namespace rpt = qchem::report;
    itsPiloted=true;
    const bool   stated=itsSetup->meshStated;
    const ClientConsumes c=Consumes();
    const bool   onE=c==ClientConsumes::Energy;
    const double target=itsSetup->opts.accuracy.Target(c, itsStructure->GetNumAtoms());
    const qcMesh::BeckeLadderCalibration cal{.angRot=0.0};   // the ladder's rungs unrotated, as the entry rung is (the rotation is a cubic-site calibration)
    const qcMesh::MeshParams entry=itsSetup->mesh;
    const std::string tag="[mesh pilot] ";
    const std::streamsize prec0=std::cout.precision();
    report::FormatGuard keep(std::cout);
    rpt::Section pilot("xcMeshPilot");
    rpt::Set("measure", onE ? "E" : "H"); rpt::Set("target", target); rpt::Set("stated", stated); rpt::Set("entry", MeshName(entry));
    const bool climbHartree=HartreeConvergesWithMesh(itsSetup->opts.coulomb);
    if (!climbHartree)
    {
        std::cout<<tag<<"WARNING: the δ-ρ DIRECT-SUM Hartree is not climbed: its error (about 1e-2) does not fall with the mesh; the climb sizes the "
                 <<"other mesh-borne terms (XC).  Leave the solver unstated (radial Poisson) for a Hartree that reaches the target"<<std::endl;
        policy::Report("DirectSumHartreeNotClimbed", "the mesh climb leaves the δ-ρ direct-sum Hartree term out of its distance", policy::Verdict::Warn,
                       "the direct sum's J error sits near 1e-2 on every rung (no better at 80 radial points than at 50), so no rung meets the target "
                       +std::to_string(target)+" on it, and chasing it costs O(N_points^2) per rung; the solver is the run's stated choice, so the "
                       "run keeps it and its Hartree is the solver's, not the mesh's");
    }
    rpt::Set("climbsHartree", climbHartree);
    std::string endWhy;
    for (int round=0; ; round++)
    {
        rpt::Row row("rounds");
        // Each rung sampled ONCE at this D: every term's energy and every (term, channel) block, keyed by the rung's spelling.
        struct Probe { std::vector<double> E; std::vector<bool> climbed; std::vector<std::vector<hmat_t<double>>> H; };   // H[term][channel]
        const qchem::ChargeDensity::rDM_CD* d=itsDensity.get();
        std::map<std::string, Probe> probes;
        auto sample=[&](const qchem::Hamiltonian::rHamiltonian& h)
        {
            Probe p;
            for (const auto& te : h.MakeTermEnergies(d)) { p.E.push_back(te.E); p.climbed.push_back(climbHartree || te.entries!="Eee"); }
            p.H.resize(p.E.size());
            for (Spin s : SpinIrreps(itsSetup->opts.spin))
                for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
                {
                    auto blocks=h.MakeTermMatrices(ibs, s, d);
                    if (blocks.size()!=p.E.size()) throw std::logic_error("MoleculeCalculation::PilotMesh: a block's term list differs from the energies'");
                    for (size_t t=0;t<blocks.size();t++) p.H[t].push_back(std::move(blocks[t]));
                }
            return p;
        };
        itsHam->RefreshForDensity(itsBasis, d);
        probes.emplace(MeshName(itsSetup->mesh), sample(*itsHam));
        auto probe=[&](const qcMesh::MeshParams& m) -> const Probe&
        {
            if (auto it=probes.find(MeshName(m)); it!=probes.end()) return it->second;
            // FENCED: a reference's own mesh and fit reports go to a nested run, so this run's `hamiltonian` keeps describing what it USED.
            struct Fence { explicit Fence(const std::string& n) {rpt::Begin(n);} ~Fence() {rpt::End();} } fence("mesh pilot reference "+MeshName(m));
            const std::unique_ptr<qchem::Hamiltonian::rHamiltonian> h(MakeHamiltonianOn(m));
            h->RefreshForDensity(itsBasis, d);
            ++itsPilotCost.referencesBuilt;
            return probes.emplace(MeshName(m), sample(*h)).first->second;   // std::map: references stay valid across later inserts
        };
        const qcMesh::BeckeRungDistance distance=[&](const qcMesh::MeshParams& a, const qcMesh::MeshParams& b)
        {
            const Probe& pa=probe(a);
            const Probe& pb=probe(b);
            if (pa.E.size()!=pb.E.size()) throw std::logic_error("MoleculeCalculation::PilotMesh: two meshes hand back different term lists");
            double worst=0.0;
            for (size_t t=0;t<pa.E.size();t++)
            {
                if (!pa.climbed[t]) continue;   // a Hartree the mesh cannot size (HartreeConvergesWithMesh)
                if (onE) { worst=std::max(worst, std::fabs(pa.E[t]-pb.E[t])); continue; }
                for (size_t k=0;k<pa.H[t].size();k++)
                {
                    const auto& x=pa.H[t][k]; const auto& y=pb.H[t][k];
                    for (size_t r=0;r<x.rows();r++) for (size_t q=r;q<x.columns();q++) worst=std::max(worst, std::fabs(x(r,q)-y(r,q)));
                }
            }
            return worst;
        };
        const int cap = (stated || round>0) ? 1 : 0;   // a stated mesh is only measured; a confirming round raises at most once
        const qcMesh::BeckeClimbResult climb=qcMesh::ClimbBeckeLadder(itsSetup->mesh, target, cal, false, distance, cap);
        itsPilotCost.rounds=round+1;
        rpt::json steps=rpt::json::array();
        bool over=false;
        for (const qcMesh::BeckeClimbStep& step : climb.steps)
        {
            over = over || step.over;
            std::cout<<tag<<"round "<<round<<"  "<<std::setw(8)<<std::left<<step.axis<<std::right<<" delivered "<<std::scientific<<std::setprecision(3)
                     <<step.delivered<<" (tail estimate "<<step.tailEstimate<<") on "<<MeshName(step.rung)<<" against "<<MeshName(step.reference)<<"  "
                     <<(!step.over ? "meets " : stated ? "OVER, and the mesh is STATED: honoured, " : "OVER ")<<(onE ? "|ΔE| " : "max|ΔH| ")<<"target "
                     <<target<<std::endl<<std::defaultfloat<<std::setprecision(prec0);
            steps.push_back({{"axis",step.axis},{"delivered",step.delivered},{"tailEstimate",step.tailEstimate},{"mesh",MeshName(step.rung)},
                             {"reference",MeshName(step.reference)},{"over",step.over}});
        }
        const char* end = climb.end==qcMesh::BeckeClimbResult::End::Settled ? "settled"
                        : climb.end==qcMesh::BeckeClimbResult::End::AxisEnd ? "axisEnd" : "raiseCap";
        rpt::Set("round", round); rpt::Set("from", MeshName(itsSetup->mesh)); rpt::Set("to", MeshName(climb.mesh));
        rpt::Set("raises", climb.raises); rpt::Set("end", end); rpt::Set("steps", steps);
        if (climb.end==qcMesh::BeckeClimbResult::End::AxisEnd) endWhy=climb.why;
        if (stated)
        {
            itsPilotCost.statedOver=over;
            itsPilotCost.settled=!over && climb.end==qcMesh::BeckeClimbResult::End::Settled;
            break;
        }
        std::cout<<tag<<"round "<<round<<"  "<<(round==0 ? "climb" : "confirming audit")<<" at this D: "<<MeshName(itsSetup->mesh);
        if (climb.raises>0) std::cout<<" -> "<<MeshName(climb.mesh)<<" in "<<climb.raises<<" raise(s)";
        std::cout<<"; "<<(climb.end==qcMesh::BeckeClimbResult::End::Settled ? std::string("both axes meet the target")
                       : climb.end==qcMesh::BeckeClimbResult::End::AxisEnd ? climb.why : std::string("one raise, audited at the next round's D"))<<std::endl;
        if (climb.raises==0) { itsPilotCost.settled = climb.end==qcMesh::BeckeClimbResult::End::Settled; break; }
        if (round+1>=kPilotRounds) { itsPilotCost.settled=false; endWhy="the pilot's round cap ("+std::to_string(kPilotRounds)+")"; break; }
        // RAISED: rebuild the point on the settled rung and re-converge ONCE from this D (the carry: same basis, block for block).
        auto next=std::make_shared<MoleculeSetup>(*itsSetup);
        next->mesh=climb.mesh;
        itsSetup=next;
        itsCarry=std::make_unique<const Carry>(MakeCarry());
        BuildGraph();
        AttachObserver();
        (void)itsScf->Iterate(params);
        itsPilotCost.reconvergeIterations+=itsScf->GetIterationCount();
        RebuildSampling();
        if (!itsScf->Converged()) { itsPilotCost.settled=false; endWhy="the re-convergence on "+MeshName(itsSetup->mesh)+" did not converge"; break; }
    }
    // FREEZE (§2 rule 4): a session steps every later point on this mesh.  A stated mesh stays stated (it was never the pilot's to move).
    if (!stated)
    {
        auto sized=std::make_shared<MoleculeSetup>(*itsSetup);
        sized->meshResolved=true;
        itsSetup=sized;
    }
    rpt::Set("rounds", itsPilotCost.rounds); rpt::Set("referencesBuilt", itsPilotCost.referencesBuilt);
    rpt::Set("reconvergeIterations", (long)itsPilotCost.reconvergeIterations); rpt::Set("settled", itsPilotCost.settled);
    rpt::Set("statedOver", itsPilotCost.statedOver); rpt::Set("mesh", MeshName(itsSetup->mesh));
    const std::string cost=std::to_string(itsPilotCost.rounds)+" round(s), "+std::to_string(itsPilotCost.referencesBuilt)+" reference Hamiltonian(s), "
                          +std::to_string(itsPilotCost.reconvergeIterations)+" re-convergence iteration(s)";
    if (itsPilotCost.statedOver)
    {
        std::cout<<tag<<"WARNING: the STATED XC mesh "<<MeshName(itsSetup->mesh)<<" misses the target "<<target<<"; honoured (leave the mesh unstated to let the run size it)"<<std::endl;
        policy::Report("StatedXCMeshMissesEps", "the stated Becke XC mesh "+MeshName(itsSetup->mesh)+" is kept", policy::Verdict::Warn,
                       "a stated mesh is honoured, never raised (pin 31 rule 3), and its delivered error is over the target "+std::to_string(target)
                       +" on "+(onE ? "|ΔE|" : "max|ΔH_ij|"));
    }
    else if (!stated && itsPilotCost.settled)
        policy::Report("MolecularXCMeshSettled", "the Becke XC mesh "+MeshName(itsSetup->mesh)+" (entered at "+MeshName(entry)+")", policy::Verdict::Choose,
                       "climbed at the converged density until both axes meet the target "+std::to_string(target)+" on "+(onE ? "|ΔE|" : "max|ΔH_ij|")+"; "+cost);
    else if (!stated)
    {
        std::cout<<tag<<"WARNING: the XC mesh is NOT SETTLED at "<<MeshName(itsSetup->mesh)<<": "<<endWhy<<std::endl;
        policy::Report("MolecularXCMeshNotSettled", "the Becke XC mesh stops at "+MeshName(itsSetup->mesh), policy::Verdict::Warn,
                       endWhy+"; the delivered error may be over the target "+std::to_string(target)+"; "+cost);
    }
}

//====================================================================================================
//  THE SAVED STATE (V-FACADE-STATE): the molecule's hooks and its restart
//====================================================================================================
namespace
{
std::string MoleculeModelName(const MoleculeSetup& su)
{
    using qchem::Hamiltonian::Model;
    if (su.pp) return "LDA";   // the PP front door is LSDA (Resolve refuses any other model)
    switch (su.opts.model)
    {
        case Model::HF:     return "HF";     case Model::E1:  return "E1";
        case Model::DHF:    return "DHF";    case Model::DE1: return "DE1";
        case Model::Xalpha: return "Xalpha"; case Model::LDA: return "LDA";
    }
    return "?";
}
} //anonymous namespace

Calculation::StateIdentity MoleculeCalculation::MakeStateIdentity() const
{
    return {.Nelec=itsSetup->Ne, .spin=itsSetup->opts.spin, .multiplicity=itsSetup->opts.multiplicity,
            .functional=MoleculeModelName(*itsSetup), .kT=itsLastkT};
}

std::vector<StateBlock> MoleculeCalculation::StateBlocks(bool withOrbitals) const
{
    if (!itsBasis) throw std::logic_error("qchem::MoleculeCalculation: no basis yet (a saved state needs the graph)");
    const auto* wf = withOrbitals && itsScf ? itsScf->GetWaveFunction() : nullptr;
    if (withOrbitals && !wf) throw std::logic_error("qchem::MoleculeCalculation: no wave function yet (a saved state needs a Converge)");
    std::vector<StateBlock> out;   // the carry's order: spin irreps outer, basis blocks inner
    for (Spin s : SpinIrreps(itsSetup->opts.spin))
        for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            StateBlock b;
            b.real=ibs;
            b.irrep=ibs->GetIrrep(s);
            if (wf) b.orbitals=wf->GetOrbitals(b.irrep);
            out.push_back(std::move(b));
        }
    return out;
}

StateSummary MoleculeCalculation::MakeStateSummary() const
{
    if (!itsScf || !itsScf->GetWaveFunction()) throw std::logic_error("qchem::MoleculeCalculation: no wave function yet");
    StateSummary s{.label="molecule", .converged=itsScf->Converged(), .iterations=itsScf->GetIterationCount(),
                   .energy=itsScf->GetEnergy().GetTotalEnergy(), .charge=itsDensity ? itsDensity->GetTotalCharge() : 0.0};
    s.recipe="model "+MoleculeModelName(*itsSetup)+", basis "+itsSetup->opts.basis+", accelerator "+itsSetup->acc.type;
    if (itsLastRow)
    {
        s.commutator   =itsLastRow->ValueOf("[F,D]");
        s.densityChange=itsLastRow->ValueOf(itsHam ? itsHam->DeltaRhoLabel() : std::string("Δρ"));
    }
    return s;
}

void MoleculeCalculation::RecordState(StateRecord& r, double kT) const
{
    Calculation::RecordState(r, kT);   // the structure, electron and block keys, and the blocks
    const CalcOptions& o=itsSetup->opts;
    StateFingerprint& fp=r.fp;
    // ---- REFUSE: the saved D is a matrix over THIS basis's functions ----
    fp.refuseStr["basis.family"] ="gaussian";
    fp.refuseStr["basis.name"]   =o.basis;
    if (!o.basisFile.empty())   // a stated file: its CONTENTS identify the basis (an edited file refuses; the path may move)
        fp.refuseStr["basis.fileChecksum"]=BasisSet::Gaussian::BasisFileChecksum(o.basisFile);
    fp.refuseStr["basis.engine"] =o.engine==Engine::LibCint ? "libcint" : "mnd";
    fp.refuseStr["basis.angular"]=o.angular==Angular::Spherical ? "spherical" : "cartesian";
    fp.refuseStr["basis.blocks"] =o.symmetry ? "SALC" : "AO";
    fp.refuseStr["pp.set"]       =itsSetup->pp ? itsSetup->pp->Name() : std::string();
    {
        std::string ov;
        for (const auto& [el,q] : o.ppOverride) ov += (ov.empty() ? "" : ",") + el + ":" + std::to_string(q);
        fp.refuseStr["pp.override"]=ov;
    }
    // The shells, from the AO basis the run built: an AO run's own blocks, a SALC run's unsymmetrized basis (a SALC block is a
    // combination of its shells and cannot report them itself).  So a .bsd edited under the same name refuses, SALC or not (lane A,
    // 2026-11-10).  A block that cannot report its shells is identified by name alone, and says so.
    {
        const BasisSet::Real_BS& ao = itsAoBasis ? *itsAoBasis : *itsBasis;
        if (o.symmetry && !itsAoBasis) throw std::logic_error("MoleculeCalculation::RecordState: a SALC run without its AO basis");
        std::vector<const BasisSet::Orbital_1E_IBS<double>*> blocks;
        bool all=true;
        for (auto ibs : ao.Iterate<BasisSet::Real_OIBS>())
        {
            blocks.push_back(ibs);
            all = all && dynamic_cast<const BasisSet::AoShellSource*>(ibs)!=nullptr;
        }
        if (all) AddShells(fp, blocks);
        else     fp.refuseStr["basis.shells"]="unreported";
    }
    // ---- WARM: what moves the converged state without changing what D is a matrix over ----
    if (o.model==Model::Xalpha) fp.warmNum["xc.xalpha"]={o.xalpha};
    if (HasXCMesh())
    {
        fp.warmNum["grids.xcMeshRadial"]={double(itsSetup->mesh.nRadial)};
        fp.warmNum["grids.xcMeshDegree"]={double(itsSetup->mesh.angularDegree)};
        fp.warmNum["accuracy.eps"]={o.accuracy.eps};
        fp.warmNum["xc.fit"]={double(int(o.vxcFit))};
    }
    fp.warmNum["coulomb.route"]={double(int(o.coulomb.fit)), double(int(o.coulomb.solver)), double(o.coulomb.lMax), double(o.coulomb.nRadial)};
}

const std::optional<RestartSummary>& MoleculeCalculation::RestartedFrom() const {return itsRestartedFrom;}

Outcome<std::unique_ptr<MoleculeCalculation>, RestartRefusal>
MoleculeCalculation::Restart(const std::string& path, std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf, Observer onIteration)
{
    auto R=PrepareRestart(path, std::move(setup), scf, std::move(onIteration));
    if (R) (void)R.Value()->Converge(scf);   // the restored stage, on the graph PrepareRestart built
    return R;
}

Outcome<std::unique_ptr<MoleculeCalculation>, RestartRefusal>
MoleculeCalculation::PrepareRestart(const std::string& path, std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf,
                                    Observer onIteration)
{
    using O=Outcome<std::unique_ptr<MoleculeCalculation>, RestartRefusal>;
    namespace rpt = qchem::report;
    const std::string tag="[molecule restart] ";
    auto saved=ReadState(path);
    if (!saved)
    {
        std::cout<<tag<<"REFUSED ("<<path<<"): "<<saved.Error().details<<std::endl;
        return O::Fail(saved.Error());
    }
    // A STATE SIZED FOR THIS TARGET (the solid's raster adoption, for the molecule's XC mesh): a setup whose mesh is neither stated
    // nor sized adopts the saved state's when the state was made for the same eps -- that mesh is what a pilot settled for this eps,
    // so an exact resume stays exact (a fresh pilot would start at the entry rung, a warm difference, and climb again).
    if (!setup->meshStated && !setup->meshResolved && saved->kind=="molecule")
    {
        const auto& w=saved->fingerprint.warmNum;
        auto one=[&](const char* k)->std::optional<double>
        { auto it=w.find(k); return it!=w.end() && it->second.size()==1 ? std::optional<double>(it->second[0]) : std::nullopt; };
        const auto nR=one("grids.xcMeshRadial"), L=one("grids.xcMeshDegree"), eps=one("accuracy.eps");
        if (nR && L && eps && *eps==setup->opts.accuracy.eps)
        {
            auto adopted=std::make_shared<MoleculeSetup>(*setup);
            adopted->mesh.nRadial=int(*nR);
            adopted->mesh.angularDegree=int(*L);
            adopted->meshResolved=true;
            std::ostringstream why;
            why<<"the saved state's XC mesh ("<<int(*nR)<<" radial, degree "<<int(*L)<<") was sized for this accuracy (eps "<<*eps<<")";
            qchem::policy::Report("RestartAdoptsMesh", "the restart integrates on the saved state's XC mesh, no pilot",
                                  qchem::policy::Verdict::Choose, why.str());
            std::cout<<tag<<"adopting "<<why.str()<<std::endl;
            setup=std::move(adopted);
        }
    }
    const Structure& ref=*setup->reference;
    std::unique_ptr<MoleculeCalculation> c(new MoleculeCalculation(BuildOnly{}, setup, ref, nullptr, std::move(onIteration)));
    {
        rpt::Scope basis("basis");
        c->itsBasis=BuildBasis(c->itsSetup->opts, c->itsStructure, c->itsAoBasis);   // BuildGraph reuses it
    }
    SavedState state=saved.TakeValue();
    auto judged=c->JudgeRestart(state, scf.SmearingkT, tag);
    if (!judged) return O::Fail(judged.Error());
    auto D=CheckedBlockDensities(state, c->StateBlocks(false));
    if (!D)
    {
        std::cout<<tag<<"REFUSED ("<<path<<"): "<<D.Error().details<<std::endl;
        return O::Fail(D.Error());
    }
    Carry carry;
    for (const mat_t<dcmplx>& M : *D)
    {
        const size_t n=M.rows();
        hmat_t<double> h(n);
        for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) h(a,b)=M(a,b).real();
        carry.D.push_back(std::move(h));
    }
    c->itsCarry=std::make_unique<const Carry>(std::move(carry));
    c->itsCarryExact=judged->exact;
    // At the answer: a Ladder starts on its GDM rung (kT = 0 only: a smeared stage restored from an integer state would let GDM
    // engage and never apply the smearing -- the solid's rule, SolidCalculation::PrepareRestart).
    c->itsStartOnGDM=judged->exact && scf.SmearingkT==0.0 && c->itsSetup->accelerator==SCFAccelerators::Type::Ladder;
    c->itsRestartedFrom=*judged;
    c->BuildGraph();
    c->itsPrepared=true;
    return O::Ok(std::move(c));
}

//====================================================================================================
//  THE MOLECULAR SESSION (doc/CalcNetworkReview.md §7, step 8)
//====================================================================================================
MoleculeSession::MoleculeSession(std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf, MoleculeCalculation::Observer onIteration)
    : itsSetup(std::move(setup)), itsSCF(scf), itsObserver(std::move(onIteration)) {}
MoleculeSession::~MoleculeSession() = default;

PESGeometry MoleculeSession::Reference() const
{
    PESGeometry g;   // no cell
    itsSetup->reference->ForEachSite([&](int Z, const rvec3_t& r, bool){ g.R.push_back(r); g.Z.push_back(Z); });
    return g;
}

Outcome<PESPoint, PESFailure> MoleculeSession::ConvergeAt(const std::optional<rmat3d_t>& cell, const std::vector<rvec3_t>& R,
                                                         const PESRequest& want)
{
    using O=Outcome<PESPoint, PESFailure>;
    using Why=PESFailure::Why;
    if (cell) throw std::invalid_argument("MoleculeSession::ConvergeAt: a molecule has no cell -- pass an empty one (Reference().cell)");
    const Structure& ref=*itsSetup->reference;
    const size_t nAtoms=ref.GetNumAtoms();
    if (R.size()!=nAtoms)
        throw std::invalid_argument("MoleculeSession::ConvergeAt: "+std::to_string(R.size())+" positions for "+std::to_string(nAtoms)+" atoms");
    if (want.forces || want.stress)
        return O::Fail({Why::NoGradient, std::string(want.stress ? "stress" : "forces")+": no term implements its geometry gradient yet "
                        "(doc/CalcNetworkReview.md §3.3, step 9) -- difference ConvergeAt energies meanwhile"});

    std::shared_ptr<Structure> at=ref.Clone();
    at->MoveSites(R);

    // THE IMPOSED GROUP must survive (Q4): the SALC blocks are built on the abelian point group, so a step that changes it changes
    // the blocks the carry and the energy surface are written in -- a NEW setup, the user's choice.  A THROW: every driver that moves
    // atoms off the group runs a lower-group setup, so an off-group R here is a wiring defect (PotentialEnergySurface.C).
    if (itsSetup->opts.symmetry)
    {
        auto abelian=[&](const Structure& st)
        {
            std::vector<qchem::Symmetry::Molecule::SymPoint> pts;
            st.ForEachSite([&](int Z, const rvec3_t& r, bool){ pts.push_back({Z, r}); });
            return qchem::Symmetry::Molecule::DetectPointGroup(pts, itsSetup->opts.symmetryTol).abelian;
        };
        if (const std::string now=abelian(*at), was=abelian(ref); now!=was)
            throw std::invalid_argument("MoleculeSession::ConvergeAt: these positions lower the imposed abelian group "+was+" to "+now
                                        +" -- step through a map given the group's operations (CartesianMap), or resolve a new setup without the symmetry");
    }

    // THE CARRY (lane A, 2026-11-09): the previous point's state, carried onto this point's basis when it is built, by the function a
    // file restart at moved atoms uses (CarryToMovedAtoms: AO coefficients, the measured best first iterate).
    std::unique_ptr<const MoleculeCalculation::Carry> carry;
    if (itsPoint)
    {
        std::vector<rvec3_t> was;
        itsPoint->StateStructure().ForEachSite([&](int, const rvec3_t& r, bool){ was.push_back(r); });
        double dmax=0.0;
        for (size_t i=0;i<R.size();i++) dmax=std::max(dmax, norm(R[i]-was[i]));
        carry=std::make_unique<const MoleculeCalculation::Carry>(itsPoint->MakeMovedCarry(dmax));
    }
    itsPoint.reset();   // one point's memory at a time
    itsPoint.reset(new MoleculeCalculation(itsSetup, *at, std::move(carry), itsSCF, itsObserver, /*geometryPoint*/true));
    auto out=itsPoint->Result();
    if (!out) return O::Fail({Why::SCF, out.Error().details, out.Error()});
    itsSetup=itsPoint->SizedSetup();   // the first point's pilot froze the mesh: every later point steps on it (§2 rule 4)
    return O::Ok(PESPoint{.energy=out->Energy(), .forces={}, .iterations=out->IterationCount()});
}

double                 MoleculeCalculation::Energy()      const {return itsScf->GetEnergy().GetTotalEnergy();}
qchem::EnergyBreakdown MoleculeCalculation::EnergyTerms() const {return itsScf->GetEnergy();}

const MoleculeCalculation::sf_t& MoleculeCalculation::Density() const
{
    assert(itsDensity);
    return *itsDensity;
}

const MoleculeCalculation::sf_t& MoleculeCalculation::HOMO() const
{
    assert(!itsOccupied.empty());
    return *itsOccupied.back().second;
}

const std::vector<MoleculeCalculation::Level>& MoleculeCalculation::GetLevels() const {return itsLevels;}

const MoleculeCalculation::sf_t& MoleculeCalculation::GetSpinDensity() const
{
    assert(itsSpinDensity);
    return *itsSpinDensity;
}

std::unique_ptr<MoleculeCalculation::sf_t> MoleculeCalculation::MakeIterateDensity() const
{
    if (!itsScf || itsScf->GetIterationCount()==0)
        throw std::logic_error("MoleculeCalculation::MakeIterateDensity: no iteration has run yet (a prepared run samples from "
                               "its first observer row on)");
    return itsScf->GetWaveFunction()->GetChargeDensity();
}

std::unique_ptr<MoleculeCalculation::sf_t> MoleculeCalculation::MakeIterateSpinDensity() const
{
    if (!itsScf || itsScf->GetIterationCount()==0)
        throw std::logic_error("MoleculeCalculation::MakeIterateSpinDensity: no iteration has run yet");
    const auto* wf=itsScf->GetWaveFunction();
    if (wf->GetSpinGroup()==SpinGroup::Polarized) return wf->GetSpinDensity();
    return std::make_unique<ZeroField>();
}

const MoleculeCalculation::sf_t& MoleculeCalculation::Orbital(size_t i) const
{
    assert(i < itsOccupied.size());
    return *itsOccupied[i].second;
}

size_t MoleculeCalculation::NumOccupied() const {return itsOccupied.size();}

const MoleculeCalculation::orbitals_t* MoleculeCalculation::Orbitals(const Irrep& qns) const
{
    return itsScf->GetWaveFunction()->GetOrbitals(qns);
}

size_t MoleculeCalculation::IterationCount() const {return itsScf->GetIterationCount();}
bool   MoleculeCalculation::IsConverged()    const {return itsScf->Converged();}

} // namespace qchem
