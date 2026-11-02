// File: Calculation/Imp/Calculation.C  Implementation of the qchem::Calculation facade.
module;
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <set>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <nlohmann/json.hpp>
module qchem.Calculation;

import qchem.BasisSet.Gaussian.Point.Factory;             // BasisSet::Gaussian::Factory
import qchem.BasisSet.IntegralsCacheReport;          // EmitIntegralsCacheReport (the run's `cache` section)
import qchem.BasisSet.Gaussian.Point.SymmetryAdaptedBasisSet; // SymmetryAdaptedBasisSet (return of SymmetryAdapt)
import qchem.BasisSet.Gaussian.PG_Cart.SymmetryAdapt;   // PG::SymmetryAdapt (the SALC builder)
import qchem.ElectronConfiguration.Molecule;        // Molecule_EC (global aufbau)
import qchem.SCFAccelerator.Factory;                // SCFAccelerators::Type, Factory
import qchem.WaveFunction;                           // WaveFunction (GetChargeDensity/GetOrbitals/GetQNs)
import qchem.Orbitals;                               // Orbital, Orbitals
import qchem.Reporting;                              // report:: run sink + the scf-section writer
import qchem.BasisSet.Orbital_1E_IBS;                // Real_OIBS -- the per-irrep blocks a carry is re-homed onto
import qchem.CompositeCD;                            // tComposite_CD (the carried density's shape)
import qchem.ChargeDensity.Factory;                  // IrrepCD_Factory
import qchem.Symmetry.Molecule.PointGroup;           // DetectPointGroup (a session step must keep the SALC group)

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
static BasisSet::Real_BS* BuildBasis(const CalcOptions& opts, const std::shared_ptr<const Structure>& st)
{
    const bool spherical = (opts.angular == Angular::Spherical);
    // SALC now supports the Cartesian PGData basis AND the in-house MnD-spherical basis (SphData;
    // doc/OldPlans/SphericalSALCPlan.md S3a/S4).  Only libcint-spherical is unwired (S3b) -- and it presents as a
    // PGData with spherical components, which the Cartesian extractor would silently misread, so reject that
    // one combination clearly here rather than let SymmetryAdapt take the wrong path.
    if (opts.symmetry && spherical && opts.engine == Engine::LibCint)
        throw std::runtime_error("qchem::Calculation: {.symmetry=true} with angular=Spherical requires "
                                 "engine=MnD (libcint-spherical SALC is not yet supported)");

    const char* engine  = (opts.engine  == Engine::LibCint) ? "libcint" : "mnd";
    const char* angular = spherical ? "spherical" : "cartesian";
    BasisSet::Real_BS* raw = BasisSet::Gaussian::Factory(
        nlohmann::json{{"basis", opts.basis}, {"engine", engine}, {"angular", angular}}, st.get());
    if (!opts.symmetry) return raw;
    std::shared_ptr<const BasisSet::Real_BS> rawShared(raw);
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

// Build the report's `scf` section json straight from the SCF knobs -- the report IS json, so this
// provider owns the keys (no mediating struct).  standard = the settled recipe surface; advanced = the
// pathological-case dials (doc/RunReportPlan.md schema).  mixer is DERIVED (Pulay > Kerker > linear);
// accelerator is the chosen engine tag.  EmitSection records it (and renders it iff a console is set).
static void EmitScfSection(const SCFParams& p, const std::string& accelerator)
{
    namespace rpt = qchem::report;
    rpt::json standard;
    standard["nMaxIter"]    = long(p.NMaxIter);
    auto gate=[](const Gate& g){return g ? rpt::json(*g) : rpt::json("off");};   // an off gate is spelled "off" (was the 1e30 sentinel)
    standard["minDrho"]     = gate(p.MinΔρPerElectron);
    standard["minDFD"]      = gate(p.MinΔFD);
    standard["minDE"]       = gate(p.MinΔE);
    standard["minFD"]       = gate(p.MinFD);
    // The preconditioner (Kerker/linear) and the history (Pulay) COMPOSE, so name both: "Ker+Pul", "Lin+Pul".
    standard["mixer"]       = std::string(p.KerkerG0 > 0.0 ? "Ker" : "Lin") + (p.PulayDepth > 0 ? "+Pul" : "");
    standard["kerkerG0"]    = p.KerkerG0;
    standard["accelerator"] = accelerator;
    standard["smearingkT"]  = p.SmearingkT;   // ortho fills when the basis/LASolver section lands

    rpt::json advanced;
    advanced["pulayDepth"]      = p.PulayDepth;
    advanced["pulayStart"]      = p.PulayStart;
    advanced["useMOM"]          = p.UseMOM;
    advanced["momStartIter"]    = p.MOMStartIter;
    advanced["momSmearPenalty"] = p.MOMSmearPenalty;
    advanced["guard"]["holePersistence"] = p.Guard.HolePersistence;
    advanced["guard"]["maxReleases"]     = p.Guard.MaxReleases;

    rpt::json scf;
    scf["standard"] = standard;
    scf["advanced"] = advanced;
    rpt::EmitSection("scf", scf);
}

Calculation::Calculation(std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf, Observer onIteration)
    : Calculation(setup, *setup->reference, nullptr, scf, std::move(onIteration))
{}

Calculation::Calculation(std::shared_ptr<const MoleculeSetup> setup, const Structure& at, std::unique_ptr<const Carry> carry,
                         const SCFParams& scf, Observer onIteration)
    : itsSetup(std::move(setup))
    // A pseudopotential run works on the VALENCE ions (Z-Zion charge) so the electron count is the valence count -- the
    // structure the basis and the EC both read.
    , itsStructure(itsSetup->PointStructure(at))
    , itsCarry(std::move(carry))
    , itsObserver(std::move(onIteration))
{
    itsEC    = MakeMoleculeEC(itsSetup->Ne, itsSetup->opts.multiplicity);
    // The basis is NOT built here: the first Converge builds it INSIDE its run report's "basis" section, so a
    // basis with something to announce does so as it is constructed (V1.14).
    (void)Converge(scf);                              // so Energy()/Density() are ready on return; Result() is its outcome
}

Calculation::~Calculation()
{
    delete itsScf;     // FIRST: the iterator references itsHam/itsAccel, which are freed after this body (R2.22)
    delete itsBasis;
    delete itsEC;
}

Outcome<Calculation::Converged, SCFFailure> Calculation::Converge(const SCFParams& params)
{
    const CalcOptions& opts=itsSetup->opts;
    nlohmann::json jsacc = {{"NProj", itsSetup->acc.nProj}, {"EMax", itsSetup->eMax},
                            {"EMin", itsSetup->acc.eMin}, {"SVTol", itsSetup->acc.svTol}};
    auto* accel = qchem::SCFAccelerators::Factory(itsSetup->accelerator, jsacc);
    const auto seed = itsSetup->seed;

    // Open a run report for the whole SCF (the RAII guard guarantees End() even if a step throws,
    // so a failed run never leaks the sink's run-context / cursor stacks).  The report only RECORDS
    // unless a client attached a console via report::SetConsole.
    namespace rpt = qchem::report;
    rpt::Begin(itsStructure->ID());
    struct RunScope { ~RunScope() { rpt::End(); } } runScope;
    EmitScfSection(params, itsSetup->acc.type);               // the "what recipe am I running" header

    // The `basis` section is ASSEMBLED across layers, and since V1.14 it spans the basis's whole LIFE in this
    // run: the FIRST Converge BUILDS the basis (a basis with something to announce does so as it is constructed --
    // nobody asks it to), the orchestrator stamps the scalars, and the WF build (inside `new SCFIter`) fills
    // basis.perIrrep + basis.removed via the cursor five layers down (MakeIrrepWFs owns the rows, the LASolver writes
    // the conditioning).  The HAMILTONIAN is built OUTSIDE it, in its own context "grids" (user ruling 2026-11-02:
    // the caller sets the context, the callee reports one level in) -- what it builds there is grids: the XC
    // quadrature, the PP terms' meshes, each filed one level in by the mesh builder, so a molecule and a solid file
    // the same paths.  So the basis context is entered twice: SILENTLY for the build (a Scope renders only the
    // scalars it added, and the build adds rows), then as the rendering Section for the rest -- which renders the
    // whole basis once, complete, on close.  A later Converge reuses the basis.
    {
        rpt::Scope basis("basis");
        if (!itsBasis) itsBasis = BuildBasis(opts, itsStructure);   // raw, or SALC-blocked if .symmetry
    }
    // The unified resolver turns the Model token into the concrete Hamiltonian; the DFT extras (mesh, orbital basis,
    // xalpha) are ignored for HF/1-e/Dirac.  A pseudopotential run takes the PP front door instead (LSDA valence
    // Hamiltonian: V_loc + KB projectors + Zion ion-ion in place of Ven).
    decltype(itsHam.get()) ham=nullptr;
    {
        rpt::Scope grids("grids");
        ham = itsSetup->pp
            ? qchem::Hamiltonian::Factory(opts.spin, itsStructure, *itsSetup->pp, opts.mesh, itsBasis, opts.vxcFit)
            : qchem::Hamiltonian::Factory(opts.model, opts.spin, itsStructure,
                                          opts.mesh, itsBasis, opts.xalpha, opts.vxcFit);
    }
    {
        rpt::Section basis("basis");
        rpt::Set("name",       opts.basis);
        rpt::Set("engine",     opts.engine  == Engine::LibCint  ? "libcint"   : "mnd");
        rpt::Set("angular",    opts.angular == Angular::Spherical ? "spherical" : "cartesian");
        rpt::Set("nFunctions", (long)itsBasis->GetNumFunctions());
        // R2.22: the iterator no longer deletes what it is handed, so this facade adopts the pair.  Order is
        // deliberate -- the PREVIOUS iterator dies first, then the previous Hamiltonian/accelerator it pointed
        // at are freed by these resets, which is exactly when `delete itsScf` used to free them.
        delete itsScf; itsScf = nullptr;
        itsHam.reset(ham);
        itsAccel.reset(accel);
        if (!itsCarry)
            itsScf = new SCFIter(itsBasis, itsEC, itsHam.get(), itsAccel.get(), seed, itsStructure.get());
        else
        {
            // THE CARRY (a MoleculeSession step): the previous point's D, block for block, onto THIS point's basis -- a
            // WorkingDensity of the moved Hamiltonian, consumed here once.
            using namespace qchem::ChargeDensity;
            auto cd=std::make_unique<tComposite_CD<double>>(itsBasis->GetReciprocalPointOps());
            size_t b=0;
            for (Spin s : SpinIrreps(opts.spin))
                for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
                {
                    if (b>=itsCarry->D.size() || itsCarry->D[b].rows()!=ibs->GetNumFunctions())
                        throw std::logic_error("qchem::Calculation: the carried density's blocks do not match this point's basis");
                    const Irrep irr=ibs->GetIrrep(s);
                    cd->Insert(std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(itsCarry->D[b++], ibs, irr)), irr);
                }
            if (b!=itsCarry->D.size())
                throw std::logic_error("qchem::Calculation: the carried density has more blocks than this point's basis");
            itsScf = new SCFIter(itsBasis, itsEC, itsHam.get(), itsAccel.get(), cd.release(), itsStructure.get(),
                                 qchem::Auto, 0.0, qchem::SCFIterator::SeedRole::WorkingDensity);
            itsCarry.reset();
        }
    }

    if (itsObserver) itsScf->SetObserver(itsObserver);
    (void)itsScf->Iterate(params);   // the outcome is the iterator's Converged(), minted below

    // Snapshot the integrals cache into the run's `cache` section (RunReportPlan step 4), inside the bracket.
    // The cache is a process-wide singleton NEVER cleared between runs, so this is cumulative-to-this-point --
    // exact for a one-run process, the running total in a multi-run one.
    qchem::BasisSet::EmitIntegralsCacheReport();

    RebuildSampling();
    return Outcome_();
}

Outcome<Calculation::Converged, SCFFailure> Calculation::Result() const {return Outcome_();}

Outcome<Calculation::Converged, SCFFailure> Calculation::Outcome_() const
{
    using O=Outcome<Converged, SCFFailure>;
    if (itsScf->Converged()) return O::Ok(Converged(this));
    SCFFailure f;   // the molecular facade's detectors: the iterator's own refusal, else a run out of iterations
    const bool refused=!itsScf->UnphysicalStop().empty();
    f.why        = refused ? SCFFailure::Why::UnphysicalEnergy : SCFFailure::Why::NotConverged;
    f.iterations = itsScf->GetIterationCount();
    f.lastEnergy = itsScf->GetEnergy().GetTotalEnergy();
    f.details    = refused ? itsScf->UnphysicalStop() : "the SCF did not converge in "+std::to_string(f.iterations)+" iterations";
    return O::Fail(std::move(f));
}

double                 Calculation::Converged::Energy()         const {return itsCalc->itsScf->GetEnergy().GetTotalEnergy();}
qchem::EnergyBreakdown Calculation::Converged::EnergyTerms()    const {return itsCalc->itsScf->GetEnergy();}
size_t                 Calculation::Converged::IterationCount() const {return itsCalc->itsScf->GetIterationCount();}
const Calculation::sf_t& Calculation::Converged::Density()      const {return itsCalc->Density();}

// The linear response of the converged state: Reference (R0, orbital basis) + frame (AO <-> MO) + the
// Hamiltonian's kernel + the dipole probe, solved together by GMRES (qcResponse).  The facade is the one place
// that owns all four ingredients, which is why the entry point is here.
Outcome<rmat_t,Response::ResponseFailure> Calculation::StaticPolarizability(const qcMesh::MeshParams& dipoleMesh,
                                                                            const KrylovParams& krylov) const
{
    const auto* wf=itsScf ? itsScf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("Calculation::StaticPolarizability: no converged wave function");
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

void Calculation::RebuildSampling()
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

void Calculation::OnIteration(Observer obs) {itsObserver = std::move(obs);}

Calculation::Carry Calculation::MakeCarry() const
{
    const auto* wf = itsScf ? itsScf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("qchem::Calculation::MakeCarry: no wave function yet");
    Carry c;
    for (Spin s : SpinIrreps(itsSetup->opts.spin))
        for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            const size_t n=ibs->GetNumFunctions();
            const auto* os=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(wf->GetOrbitals(ibs->GetIrrep(s)));
            if (!os) throw std::logic_error("qchem::Calculation::MakeCarry: a block's orbitals are not real");
            hmat_t<double> D(n);
            for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) D(a,b)=0.0;
            for (const auto* o : os->Iterate<qchem::Orbitals::TOrbital<double>>())
            {
                const double f=o->GetOccupation();
                if (f==0.0) continue;
                const auto& v=o->GetCoeff();
                for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) D(a,b)+=f*v[a]*v[b];
            }
            c.D.push_back(std::move(D));
        }
    return c;
}

//====================================================================================================
//  THE MOLECULAR SESSION (doc/CalcNetworkReview.md §7, step 8)
//====================================================================================================
MoleculeSession::MoleculeSession(std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf, Calculation::Observer onIteration)
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
    for (size_t a=0;a<nAtoms;a++) (*at)[a]->itsR=R[a];

    // THE IMPOSED GROUP must survive (Q4): the SALC blocks are built on the abelian point group, so a step that changes it changes
    // the blocks the carry and the energy surface are written in -- a NEW setup, the user's choice.
    if (itsSetup->opts.symmetry)
    {
        auto abelian=[&](const Structure& st)
        {
            std::vector<qchem::Symmetry::Molecule::SymPoint> pts;
            st.ForEachSite([&](int Z, const rvec3_t& r, bool){ pts.push_back({Z, r}); });
            return qchem::Symmetry::Molecule::DetectPointGroup(pts, itsSetup->opts.symmetryTol).abelian;
        };
        if (const std::string now=abelian(*at), was=abelian(ref); now!=was)
            return O::Fail({Why::SymmetryBroken, "these positions lower the imposed abelian group "+was+" to "+now
                            +" -- symmetrize the step, or resolve a new setup without the symmetry"});
    }

    std::unique_ptr<const Calculation::Carry> carry;
    if (itsPoint) carry=std::make_unique<const Calculation::Carry>(itsPoint->MakeCarry());
    itsPoint.reset();   // one point's memory at a time
    itsPoint.reset(new Calculation(itsSetup, *at, std::move(carry), itsSCF, itsObserver));
    auto out=itsPoint->Result();
    if (!out) return O::Fail({Why::SCF, out.Error().details, out.Error()});
    return O::Ok(PESPoint{.energy=out->Energy(), .forces={}, .iterations=out->IterationCount()});
}

double                 Calculation::Energy()      const {return itsScf->GetEnergy().GetTotalEnergy();}
qchem::EnergyBreakdown Calculation::EnergyTerms() const {return itsScf->GetEnergy();}

const Calculation::sf_t& Calculation::Density() const
{
    assert(itsDensity);
    return *itsDensity;
}

const Calculation::sf_t& Calculation::HOMO() const
{
    assert(!itsOccupied.empty());
    return *itsOccupied.back().second;
}

const std::vector<Calculation::Level>& Calculation::GetLevels() const {return itsLevels;}

const Calculation::sf_t& Calculation::GetSpinDensity() const
{
    assert(itsSpinDensity);
    return *itsSpinDensity;
}

const Calculation::sf_t& Calculation::Orbital(size_t i) const
{
    assert(i < itsOccupied.size());
    return *itsOccupied[i].second;
}

size_t Calculation::NumOccupied() const {return itsOccupied.size();}

const Calculation::orbitals_t* Calculation::Orbitals(const Irrep& qns) const
{
    return itsScf->GetWaveFunction()->GetOrbitals(qns);
}

size_t Calculation::IterationCount() const {return itsScf->GetIterationCount();}
bool   Calculation::IsConverged()    const {return itsScf->Converged();}

} // namespace qchem
