// File: pybind/qchem_bridge.cpp
//
// Module unit: the global-module-fragment (#includes) de-dups against the imported qchem modules' GMFs, so
// textually including STL alongside `import qchem.*` does not trigger std-header ODR errors.  The public
// surface is the flat extern "C" API in qcb_api.h (C linkage -> links from nanobind).
//
// The bridge holds NO recipe of its own.  A molecule: it fills CalcOptions / AcceleratorOptions / SCFParams, calls
// MoleculeSetup::Resolve, and owns one qchem::MoleculeCalculation.  A crystal: it parses the caller's INPUT DECK and calls
// deck::Build, owning the SolidCalculation that returns.  Everything it answers comes from the facade.
module;
#include <algorithm>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <cassert>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "qcb_api.h"

export module qchem.bridge;

import qchem.Types;          // rvec3_t, vec3_t
import qchem.Structure;      // Molecule, Site
import qchem.SCFParams;      // SCFParams
import qchem.SCFIterator;    // IterationRow, IterationCell
import qchem.MoleculeCalculation;    // MoleculeCalculation, MoleculeSetup, CalcOptions, AcceleratorOptions
import qchem.SolidCalculation;  // SolidCalculation (what deck::Build returns)
import qchem.Deck;           // deck::RunSpec, FromJson, Build, ScheduleToRun
import qchem.Symmetry.Molecule.PointGroup;      // DetectPointGroup
import qchem.Symmetry.Lattice_3D.SpaceGroup;    // SpaceGroup::Detect
import qchem.Reporting;                          // report::GlobalReport (the run record the library keeps)
import qchem.CrystalFiles;                       // CIF / POSCAR / extended xyz <-> one materials.json entry
import qchem.StructureData;                      // the materials layers (shipped, user, workspace) and CellFromEntry

// The opaque handle: the facade plus the request's geometry (for the sampling box).  extern "C++" attaches it to the
// GLOBAL module, where qcb_api.h (in the GMF) forward-declares it -- a purview definition would be a different entity.
extern "C++" struct qcb_calc
{
    std::unique_ptr<qchem::MoleculeCalculation> calc;
    qchem::SCFParams                    params;      // as asked at qcb_make; a re-run changes only NMaxIter
    std::vector<double>                 positions;   // flat 3N bohr, as asked
    std::string                         id;          // the structure's ID: the name of the report runs the bridge opens
};

//! A crystal: the deck as read, and the calculation deck::Build made of it.
extern "C++" struct qcb_solid
{
    qchem::deck::RunSpec                       spec;
    std::unique_ptr<qchem::SolidCalculation>   calc;
    bool                                       attempted = false;   // Result() before any Converge has nothing to say
};

namespace {

using qchem::MoleculeCalculation;
using qchem::SolidCalculation;

//! The HOST's C locale, put back when an entry point returns.  The lib's console tables (submodules/tabulate) bracket each
//! border with std::locale::global(en_US.UTF-8) ... std::locale::global(<the saved C++ global>), and the C++ global of a
//! process that never set one is "C": so every SCF left Python's process in LC_ALL=C (measured 2026-10-31 -- a subprocess
//! read with text=True then decoded as ASCII).  Every extern "C" entry that can print holds one.
class HostLocale
{
public:
    HostLocale() : itsSaved(std::setlocale(LC_ALL, nullptr)) {}
    ~HostLocale() { std::setlocale(LC_ALL, itsSaved.c_str()); }
    HostLocale(const HostLocale&)            = delete;
    HostLocale& operator=(const HostLocale&) = delete;
private:
    std::string itsSaved;
};

void write_error(const char* what, char* err, int errlen)
{
    if (!err || errlen <= 0) return;
    std::strncpy(err, what, size_t(errlen) - 1);
    err[errlen - 1] = '\0';
}

qchem::Model parse_model(const std::string& m)
{
    using qchem::Model;
    if (m == "HF")     return Model::HF;
    if (m == "LDA")    return Model::LDA;
    if (m == "Xalpha") return Model::Xalpha;
    throw std::invalid_argument("qchem_py: unknown model '" + m + "' (HF | LDA | Xalpha)");
}

qchem::SCFParams scf_params(const qcb_request& req)
{
    qchem::SCFParams p;
    p.NMaxIter = size_t(req.maxIter);
    if (req.relaxStart > 0.0) p.StartingRelaxRo = req.relaxStart;
    if (req.relaxMax   > 0.0) p.MaxRelaxRo      = req.relaxMax;
    if (req.relaxMin   > 0.0) p.MinRelaxRo      = req.relaxMin;
    const auto gate = [](qchem::Gate& g, double v) { if (v > 0.0) g = v; else if (v == 0.0) g = std::nullopt; };   // < 0: the default
    gate(p.MinFD, req.gateFD); gate(p.MinΔρPerElectron, req.gateDrho); gate(p.MinΔE, req.gateDE); gate(p.MinΔFD, req.gateDFD);
    p.UseMOM = req.useMOM != 0;
    if (req.momStart   > 0)   p.MOMStartIter    = req.momStart;
    using M = qchem::SCFParams::Mixing;
    const std::string mixer = req.mixer ? req.mixer : "";
    p.Mixer = mixer.empty() || mixer == "Auto" ? M::Auto : mixer == "Linear" ? M::Linear : mixer == "Kerker" ? M::Kerker
            : mixer == "Pulay" ? M::Pulay
            : throw std::invalid_argument("qchem_py: unknown mixer \"" + mixer + "\" (Auto | Linear | Kerker | Pulay)");
    p.SmearingkT = req.smearingkT;
    return p;
}

//! Thrown through the SCF loop when the sink asks to stop (qcb_row_fn's nonzero return); caught at the C boundary.
struct StoppedByUser : std::runtime_error { StoppedByUser() : std::runtime_error("stopped by the user") {} };

//! The facade's observer, forwarding each iteration's row to the C sink: every cell, generically.  One row type serves
//! both facades (MoleculeCalculation::Observer and SolidCalculation::Observer take the same IterationRow).
std::function<void(const qchem::SCFIterator::IterationRow&)> forward_to(qcb_row_fn step, qcb_sink* sink)
{
    return [step, sink](const qchem::SCFIterator::IterationRow& r)
    {
        using Scale = qchem::SCFIterator::IterationCell::Scale;
        std::vector<qcb_cell> cells;
        cells.reserve(r.cells.size());
        for (const auto& c : r.cells)
            cells.push_back({c.title.c_str(), c.value, c.text.c_str(), c.tolerance.value_or(std::nan("")),
                             c.scale == Scale::Log ? 1 : 0});
        const qcb_row row{int(r.iteration), r.eb.GetTotalEnergy(), int(cells.size()), cells.data()};
        if (step(sink, &row)) throw StoppedByUser{};
    };
}

std::string why_name(qchem::SCFFailure::Why w)
{
    using Why = qchem::SCFFailure::Why;
    switch (w)
    {
        case Why::NotConverged: return "not converged";
        case Why::ChargeSlosh:  return "charge slosh (the Hartree term ran away)";
        case Why::OrderLost:    return "magnetic order lost";
        case Why::UnphysicalEnergy: return "unphysical first iterate (electronic energy > 0)";
        case Why::ResourceRefused:  return "refused on resources (the run would not fit)";
    }
    return "unknown SCF failure";
}

//! The SCFFailure as one line for Python: why, after how many iterations, and the facade's details (which gates missed, and by
//! how much, are the run's own last row: the GUI shows it).
std::string describe(const qchem::SCFFailure& f)
{
    std::string s = why_name(f.why) + " after " + std::to_string(f.iterations) + " iterations";
    if (!f.details.empty()) s += ": " + f.details;
    return s;
}

void bbox(const std::vector<double>& pos, double pad, rvec3_t& lo, rvec3_t& hi)
{
    double mn[3] = { 1e30,  1e30,  1e30}, mx[3] = {-1e30, -1e30, -1e30};
    for (size_t a = 0; a < pos.size()/3; ++a)
        for (int d = 0; d < 3; ++d)
        {
            mn[d] = std::min(mn[d], pos[3*a+d]);
            mx[d] = std::max(mx[d], pos[3*a+d]);
        }
    lo = rvec3_t(mn[0]-pad, mn[1]-pad, mn[2]-pad);
    hi = rvec3_t(mx[0]+pad, mx[1]+pad, mx[2]+pad);
}

//! Visit the n^3 grid points of the padded box (C-order, k fastest); write origin/spacing.
template <class F> void for_grid(const qcb_calc* h, int n, double pad, double* origin, double* spacing, F&& f)
{
    rvec3_t lo, hi; bbox(h->positions, pad, lo, hi);
    const rvec3_t sp((hi.x-lo.x)/(n-1), (hi.y-lo.y)/(n-1), (hi.z-lo.z)/(n-1));
    origin[0]=lo.x;  origin[1]=lo.y;  origin[2]=lo.z;
    spacing[0]=sp.x; spacing[1]=sp.y; spacing[2]=sp.z;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                f(rvec3_t(lo.x+sp.x*i, lo.y+sp.y*j, lo.z+sp.z*k));
}

void sample(const qcb_calc* h, const MoleculeCalculation::sf_t& f, int n, double pad, double* out, double* origin, double* spacing)
{
    size_t idx = 0;
    for_grid(h, n, pad, origin, spacing, [&](const rvec3_t& r){ out[idx++] = f(r); });
}

//! What a molecular request resolves to: the molecule and the facade's options, as qcb_make and qcb_restart both build them.
struct MoleculeRequest
{
    qchem::Molecule           mol;
    qchem::CalcOptions        opts;
    qchem::AcceleratorOptions acc;
};

MoleculeRequest molecule_request(const qcb_request* req)
{
    MoleculeRequest r;
    for (int a = 0; a < req->nat; ++a)
        r.mol.Insert(qchem::Site{req->Z[a], rvec3_t(req->pos3[3*a], req->pos3[3*a+1], req->pos3[3*a+2])});
    r.opts.basis        = req->basis;
    r.opts.model        = parse_model(req->model);
    r.opts.pp           = req->pp ? req->pp : "";
    r.opts.symmetry     = req->symmetry != 0;
    r.opts.spin         = req->polarized ? qchem::SpinGroup::Polarized : qchem::SpinGroup::UnPolarized;
    r.opts.multiplicity = req->multiplicity;
    if (req->meshRadial  > 0) r.opts.mesh.nRadial       = req->meshRadial;
    if (req->meshAngular > 0) r.opts.mesh.angularDegree = req->meshAngular;
    r.acc.type = req->accelerator && *req->accelerator ? req->accelerator : "DIIS";
    if (req->diisDepth    > 0)   r.acc.nProj = req->diisDepth;
    if (req->diisSVTol    > 0.0) r.acc.svTol = req->diisSVTol;
    if (req->diisEngageFD > 0.0) r.acc.eMax  = req->diisEngageFD;
    if (req->ladderSwitchAt >= 0.0) r.acc.ladderSwitchAt = req->ladderSwitchAt;   // 0 = no tail trigger (stall only)
    return r;
}

const char* why_name(qchem::RestartRefusal::Why w)
{
    using W = qchem::RestartRefusal::Why;
    switch (w)
    {
        case W::Unreadable: return "unreadable";
        case W::Format:     return "format";
        case W::Kind:       return "kind";
        case W::Mismatch:   return "mismatch";
        case W::Realness:   return "realness";
    }
    return "unknown";
}

nlohmann::ordered_json refusal_json(const qchem::RestartRefusal& r)
{
    return {{"refused", {{"why", why_name(r.why)}, {"details", r.details}}}};
}

//! A facade's RestartedFrom as JSON: null for a fresh run.
nlohmann::ordered_json restarted_json(const std::optional<qchem::RestartSummary>& R)
{
    if (!R) return nullptr;
    return {{"path", R->path}, {"exact", R->exact}, {"warm", R->warm}, {"maxDisplacement", R->maxDisplacement},
            {"electrons", R->electrons}};
}

//! JudgeState as JSON: {"exact": B, "warm": [...]} or {"refused": {why, details}}.
nlohmann::ordered_json judged_json(const qchem::Outcome<std::vector<std::string>, qchem::RestartRefusal>& j)
{
    if (!j) return refusal_json(j.Error());
    return {{"exact", j.Value().empty()}, {"warm", j.Value()}};
}

int write_json(const nlohmann::ordered_json& j, char* out, int outlen)
{
    const std::string s = j.dump(-1, ' ', false, nlohmann::ordered_json::error_handler_t::replace);
    write_error(s.c_str(), out, outlen);
    return int(s.size());
}

} // anon namespace

extern "C" {

qcb_calc* qcb_make(const qcb_request* req, qcb_row_fn step, qcb_sink* sink, char* err, int errlen)
{
    return qcb_start(req, nullptr, 0, step, sink, nullptr, err, errlen);
}

qcb_calc* qcb_start(const qcb_request* req, const char* path, int buildOnly, qcb_row_fn step, qcb_sink* sink, int* refused,
                    char* err, int errlen)
{
    const HostLocale host;
    if (refused) *refused = 0;
    try
    {
        auto h = std::make_unique<qcb_calc>();
        h->positions.assign(req->pos3, req->pos3 + 3*req->nat);
        auto r = molecule_request(req);

        // The bridge is the DRIVER (user ruling D): it opens the run the facade writes into -- "<ID>", as for a Converge --
        // BEFORE the setup's checks, so a refusal there files its record (its policy rows) too.
        // A BUILD (no SCF) files "<ID> build", as a crystal's deck build does; its first qcb_converge reuses the graph.
        h->id     = r.mol.ID();
        const qchem::report::Run run(buildOnly ? h->id + " build" : h->id);
        const auto setup = qchem::MoleculeSetup::Resolve(r.mol, r.opts, r.acc);   // deep-copies mol
        // The observer rides the constructor, so the trace of ITS SCF streams from the first iteration.
        h->params = scf_params(*req);
        auto observer = step ? forward_to(step, sink) : MoleculeCalculation::Observer{};
        if (path && *path)
        {
            // A refusal is a VALUE (the facade's RestartRefusal): the caller decides on a fresh run, so it is told apart from a failure.
            auto R = buildOnly ? MoleculeCalculation::PrepareRestart(path, setup, h->params, std::move(observer))
                               : MoleculeCalculation::Restart(path, setup, h->params, std::move(observer));
            if (!R)
            {
                if (refused) *refused = 1;
                write_error(("restart refused (" + std::string(why_name(R.Error().why)) + "): " + R.Error().details).c_str(), err, errlen);
                return nullptr;
            }
            h->calc = R.TakeValue();
        }
        else if (buildOnly)
            h->calc = MoleculeCalculation::Prepare(setup, std::move(observer));
        else
            h->calc = std::make_unique<MoleculeCalculation>(setup, h->params, std::move(observer));
        // The sink lives only for this call: detach it before returning (qcb_converge attaches its own).
        h->calc->OnIteration(MoleculeCalculation::Observer{});
        return h.release();
    }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    return nullptr;
}

int qcb_save_state(const qcb_calc* h, const char* path, char* err, int errlen)
{
    const HostLocale host;
    try
    {
        // The save is its own run in the record, "<ID> state save" (its timing), so the SCF run stays the last one a reader picks.
        const qchem::report::Run run(h->id + " state save");
        h->calc->SaveState(path);
        return 1;
    }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    return 0;
}

int qcb_judge_state(const qcb_calc* h, const char* path, char* out, int outlen)
{
    const HostLocale host;
    return write_json(judged_json(h->calc->JudgeState(path)), out, outlen);
}

int qcb_restarted_from(const qcb_calc* h, char* out, int outlen)
{
    return write_json(restarted_json(h->calc->RestartedFrom()), out, outlen);
}

void qcb_free(qcb_calc* h) { delete h; }

int qcb_converge(qcb_calc* h, int maxIter, qcb_row_fn step, qcb_sink* sink, char* err, int errlen)
{
    const HostLocale host;
    try
    {
        h->calc->OnIteration(step ? forward_to(step, sink) : MoleculeCalculation::Observer{});
        qchem::SCFParams p = h->params;
        p.NMaxIter = size_t(maxIter);
        const auto out = [&]{ const qchem::report::Run run(h->id); return h->calc->Converge(p); }();
        h->calc->OnIteration(MoleculeCalculation::Observer{});      // the sink lives only for this call
        if (out) return 1;
        write_error(describe(out.Error()).c_str(), err, errlen);
        return 0;
    }
    catch (const StoppedByUser& e) { h->calc->OnIteration(MoleculeCalculation::Observer{}); write_error(e.what(), err, errlen); return 0; }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    return -1;
}

int qcb_failure(const qcb_calc* h, char* out, int outlen)
{
    const auto r = h->calc->Result();
    write_error(r ? "" : describe(r.Error()).c_str(), out, outlen);
    return r ? 0 : 1;
}

double qcb_energy    (const qcb_calc* h) { return h->calc->Energy(); }
int    qcb_converged (const qcb_calc* h) { return h->calc->IsConverged() ? 1 : 0; }
int    qcb_iterations(const qcb_calc* h) { return int(h->calc->IterationCount()); }
int    qcb_natoms    (const qcb_calc* h) { return int(h->positions.size()/3); }

void qcb_atoms(const qcb_calc* h, int* Z_out, double* xyz_out)
{
    const qchem::Structure& st = h->calc->GetStructure();
    size_t a = 0;
    for (const qchem::Site& site : st.Sites()) Z_out[a++] = site.Z;
    std::memcpy(xyz_out, h->positions.data(), h->positions.size()*sizeof(double));
}

int qcb_num_levels(const qcb_calc* h) { return int(h->calc->GetLevels().size()); }

void qcb_levels(const qcb_calc* h, double* e, double* occ, int* spin, char* label, int nlabel)
{
    const auto& levels = h->calc->GetLevels();
    for (size_t i = 0; i < levels.size(); ++i)
    {
        e[i]    = levels[i].e;
        occ[i]  = levels[i].occupation;
        spin[i] = levels[i].spin == qchem::Spin::Up ? 1 : levels[i].spin == qchem::Spin::Down ? -1 : 0;
        write_error(levels[i].label.c_str(), label + i*nlabel, nlabel);
    }
}

void qcb_density(const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing)
{
    sample(h, h->calc->Density(), n, pad, out, origin, spacing);
}

int qcb_iterate_density(const qcb_calc* h, int spin, int n, double pad, double* out, double* origin, double* spacing,
                        char* err, int errlen)
{
    const HostLocale host;
    try
    {
        const auto f = spin ? h->calc->MakeIterateSpinDensity() : h->calc->MakeIterateDensity();
        sample(h, *f, n, pad, out, origin, spacing);
        return 1;
    }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    return 0;
}

void qcb_spin_density(const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing)
{
    sample(h, h->calc->GetSpinDensity(), n, pad, out, origin, spacing);
}

int qcb_orbital(const qcb_calc* h, int index, int n, double pad, double* out, double* origin, double* spacing)
{
    const auto& levels = h->calc->GetLevels();
    if (index < 0 || size_t(index) >= levels.size()) return 0;
    sample(h, *levels[size_t(index)].phi, n, pad, out, origin, spacing);
    return 1;
}

void qcb_gradient(const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing)
{
    const auto& rho = h->calc->Density();
    size_t idx = 0;
    for_grid(h, n, pad, origin, spacing, [&](const rvec3_t& r)
    {
        const vec3_t<double> g = rho.Gradient(r);
        out[idx++] = g.x; out[idx++] = g.y; out[idx++] = g.z;
    });
}

// -- periodic runs ---------------------------------------------------------------------------------------------------

qcb_solid* qcb_solid_make(const char* deckJson, char* err, int errlen)
{
    const HostLocale host;
    try
    {
        auto h = std::make_unique<qcb_solid>();
        const auto deck = nlohmann::json::parse(deckJson);
        // The bridge is the driver: it opens the runs -- this one BEFORE the deck is read, so a refusal while reading it files
        // its record too.  Its name is RunSpec::StructureName's, read off the raw deck (the spec does not exist yet).
        const auto& st = deck.contains("structure") ? deck["structure"] : nlohmann::json();
        const std::string name = st.is_string() ? st.get<std::string>() : st.is_object() ? st.value("name", std::string("inline")) : "deck";
        const qchem::report::Run run(name + " build");
        qchem::deck::FromJson(deck, h->spec);                               // strict: an unknown key throws, naming it
        h->calc = qchem::deck::Build(h->spec);                              // the observer is attached per Converge
        return h.release();
    }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    return nullptr;
}

void qcb_solid_free(qcb_solid* h) { delete h; }

int qcb_solid_converge(qcb_solid* h, qcb_row_fn step, qcb_sink* sink, char* err, int errlen)
{
    const HostLocale host;
    try
    {
        h->calc->OnIteration(step ? forward_to(step, sink) : SolidCalculation::Observer{});
        h->attempted = true;
        const qchem::report::Run run(h->spec.StructureName());
        const auto out = h->calc->Converge(qchem::deck::ScheduleToRun(h->spec));
        h->calc->OnIteration(SolidCalculation::Observer{});      // the sink lives only for this call
        if (out) return 1;
        write_error(describe(out.Error()).c_str(), err, errlen);
        return 0;
    }
    catch (const StoppedByUser& e)  { h->calc->OnIteration(SolidCalculation::Observer{}); write_error(e.what(), err, errlen); return 0; }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    h->calc->OnIteration(SolidCalculation::Observer{});
    return -1;
}

int qcb_solid_failure(const qcb_solid* h, char* out, int outlen)
{
    if (!h->attempted) { write_error("not converged yet", out, outlen); return 1; }
    const auto r = h->calc->Result();
    write_error(r ? "" : describe(r.Error()).c_str(), out, outlen);
    return r ? 0 : 1;
}

double qcb_solid_energy(const qcb_solid* h)
{
    if (!h->attempted) return std::nan("");
    const auto r = h->calc->Result();
    return r ? r.Value().Energy() : std::nan("");
}

int qcb_solid_iterations(const qcb_solid* h) { return int(h->calc->IterationCount()); }

int qcb_solid_resolved(const qcb_solid* h, char* out, int outlen)
{
    qchem::deck::RunSpec spec = h->spec;
    (void)qchem::deck::Resolve(spec);                     // fills the derived Nelec / species, as Run's record does
    const std::string s = qchem::deck::ToJson(spec).dump(2);
    write_error(s.c_str(), out, outlen);
    return int(s.size());
}

int qcb_solid_judge_state(const qcb_solid* h, const char* path, char* out, int outlen)
{
    const HostLocale host;
    return write_json(judged_json(h->calc->JudgeState(path)), out, outlen);
}

int qcb_solid_restarted_from(const qcb_solid* h, char* out, int outlen)
{
    return write_json(restarted_json(h->calc->RestartedFrom()), out, outlen);
}

// -- the run record ---------------------------------------------------------------------------------------------------

int qcb_report(char* out, int outlen)
{
    const std::string s = qchem::report::GlobalReport().dump(-1, ' ', false, nlohmann::ordered_json::error_handler_t::replace);
    write_error(s.c_str(), out, outlen);
    return int(s.size());
}

int qcb_current_report(char* out, int outlen)
{
    const std::string s = qchem::report::Depth()>0
        ? qchem::report::CurrentRunReport().dump(-1, ' ', false, nlohmann::ordered_json::error_handler_t::replace) : std::string();
    write_error(s.c_str(), out, outlen);
    return int(s.size());
}

void qcb_clear_report() { if (qchem::report::Depth()==0) qchem::report::ClearGlobal(); }

// -- detected symmetry -----------------------------------------------------------------------------------------------

int qcb_point_group(const int* species, int nat, const double* pos3, double tol,
                    char* symbol, int symlen, char* abelian, int ablen)
{
    try
    {
        std::vector<qchem::Symmetry::Molecule::SymPoint> pts;
        for (int a = 0; a < nat; ++a)
            pts.push_back({species[a], rvec3_t(pos3[3*a], pos3[3*a+1], pos3[3*a+2])});
        const auto g = qchem::Symmetry::Molecule::DetectPointGroup(pts, tol);
        write_error(g.symbol.c_str(), symbol, symlen);
        write_error(g.abelian.c_str(), abelian, ablen);
        return g.order;
    }
    catch (const std::exception& e) { write_error(e.what(), symbol, symlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", symbol, symlen); }
    return -1;
}

int qcb_space_group(const double* L, const int* species, int nat, const double* frac3, double tol,
                    double* W9, int maxops, int* symmorphic, char* err, int errlen)
{
    try
    {
        using namespace qchem::Symmetry::Lattice_3D;
        // the cell matrix's COLUMNS are the lattice vectors; the caller's rows are
        const qchem::Matrix3D<double> A(L[0], L[3], L[6], L[1], L[4], L[7], L[2], L[5], L[8]);
        std::vector<AtomSite> basis;
        for (int a = 0; a < nat; ++a)
            basis.push_back({species[a], rvec3_t(frac3[3*a], frac3[3*a+1], frac3[3*a+2])});
        const SpaceGroup g = SpaceGroup::Detect(A, basis, tol);
        const auto& ops = g.Ops();
        for (size_t k = 0; k < ops.size() && int(k) < maxops; ++k)
            for (size_t i = 1; i <= 3; ++i)
                for (size_t j = 1; j <= 3; ++j)
                    W9[9*k + 3*(i-1) + (j-1)] = ops[k].W(i, j);
        *symmorphic = g.isSymmorphic() ? 1 : 0;
        return int(ops.size());
    }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    return -1;
}

} // extern "C"

// -- crystal files and the materials layers -------------------------------------------------------------------------------------

namespace
{
nlohmann::ordered_json LatticeRows(const std::string& entryJson)
{
    const auto A=qchem::StructureData::CellFromEntry(entryJson, "entry").GetCellMatrix();   // columns are the lattice vectors
    nlohmann::ordered_json rows=nlohmann::ordered_json::array();
    for (int i=1;i<=3;i++) rows.push_back({A(1,i), A(2,i), A(3,i)});
    return rows;
}

nlohmann::ordered_json MaterialsAnswer(const std::string& req)
{
    namespace SD=qchem::StructureData;
    namespace CF=qchem::CrystalFiles;
    using json=nlohmann::ordered_json;
    const json r=json::parse(req);
    const std::string op=r.at("op").get<std::string>();
    auto layerOf=[](const json& l)
    {
        const std::string k=l.at("layer").get<std::string>();
        if (k!="user" && k!="workspace") throw std::invalid_argument("a library layer is user or workspace, not '"+k+"'");
        return SD::LibraryFile{k=="user" ? SD::Layer::User : SD::Layer::Workspace, l.at("path").get<std::string>()};
    };
    std::vector<SD::LibraryFile> libraries;
    for (const auto& l : r.value("libraries", json::array())) libraries.push_back(layerOf(l));
    if (op=="list")
    {
        json all=json::array();
        for (const auto& e : SD::CellEntries(libraries))
            all.push_back({{"name", e.name}, {"layer", e.layer==SD::Layer::Shipped ? "shipped" : e.layer==SD::Layer::User ? "user" : "workspace"},
                           {"file", e.file}, {"entry", json::parse(e.entryJson)}, {"lattice_bohr", LatticeRows(e.entryJson)}});
        return {{"materials", all}};
    }
    if (op=="import")
    {
        const auto imp=CF::ReadFile(r.at("path").get<std::string>());
        return {{"entry", json::parse(imp.entry)}, {"name", imp.name}, {"warnings", imp.warnings}, {"lattice_bohr", LatticeRows(imp.entry)}};
    }
    if (op=="export")
    {
        const std::string f=r.at("format").get<std::string>();
        const CF::Format fmt = f=="cif" ? CF::Format::CIF : f=="poscar" ? CF::Format::POSCAR : f=="xyz" ? CF::Format::ExtendedXYZ
                             : throw std::invalid_argument("export format '"+f+"' (known: cif, poscar, xyz)");
        return {{"text", CF::Write(r.at("entry").dump(), fmt, r.at("name").get<std::string>())}};
    }
    if (op=="save")
    {
        const SD::LibraryFile into=layerOf(r.at("into"));
        SD::SaveMaterial(into, r.at("name").get<std::string>(), r.at("entry").dump(), r.value("replace", false), libraries);
        return {{"path", into.path}};
    }
    if (op=="cell") return {{"lattice_bohr", LatticeRows(r.at("entry").dump())}};
    throw std::invalid_argument("materials op '"+op+"' (known: list, import, export, save, cell)");
}
} // namespace

extern "C" int qcb_materials(const char* req, char* out, int outlen)
{
    static thread_local std::string answer;
    if (!out)
    {
        try { answer=MaterialsAnswer(req ? req : "").dump(); }
        catch (const std::exception& e) { answer=nlohmann::ordered_json{{"error", e.what()}}.dump(); }
        return int(answer.size());
    }
    write_error(answer.c_str(), out, outlen);
    return int(answer.size());
}
