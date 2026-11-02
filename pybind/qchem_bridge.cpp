// File: pybind/qchem_bridge.cpp
//
// Module unit: the global-module-fragment (#includes) de-dups against the imported qchem modules' GMFs, so
// textually including STL alongside `import qchem.*` does not trigger std-header ODR errors.  The public
// surface is the flat extern "C" API in qcb_api.h (C linkage -> links from nanobind).
//
// The bridge holds NO recipe of its own.  A molecule: it fills CalcOptions / AcceleratorOptions / SCFParams, calls
// MoleculeSetup::Resolve, and owns one qchem::Calculation.  A crystal: it parses the caller's INPUT DECK and calls
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
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "qcb_api.h"

export module qchem.bridge;

import qchem.Types;          // rvec3_t, vec3_t
import qchem.Structure;      // Molecule, Atom
import qchem.SCFParams;      // SCFParams
import qchem.SCFIterator;    // IterationRow, IterationCell
import qchem.Calculation;    // Calculation, MoleculeSetup, CalcOptions, AcceleratorOptions
import qchem.SolidCalculation;  // SolidCalculation (what deck::Build returns)
import qchem.Deck;           // deck::RunSpec, FromJson, Build, ScheduleToRun
import qchem.Symmetry.Molecule.PointGroup;      // DetectPointGroup
import qchem.Symmetry.Lattice_3D.SpaceGroup;    // SpaceGroup::Detect
import qchem.BasisSet;                           // Real_BS (Iterate over its irrep blocks)
import qchem.BasisSet.IrrepBasisSet;             // Integrals_Overlap (each block's S)
import qchem.BasisSet.Gaussian.Point.Factory;    // BasisSet::Gaussian::Factory (the facade's own basis builder)
import qchem.Blaze;                              // blazem::eigen
import qchem.Reporting;                          // report::GlobalReport (the run record the library keeps)

// The opaque handle: the facade plus the request's geometry (for the sampling box).  extern "C++" attaches it to the
// GLOBAL module, where qcb_api.h (in the GMF) forward-declares it -- a purview definition would be a different entity.
extern "C++" struct qcb_calc
{
    std::unique_ptr<qchem::Calculation> calc;
    qchem::SCFParams                    params;      // as asked at qcb_make; a re-run changes only NMaxIter
    std::vector<double>                 positions;   // flat 3N bohr, as asked
};

//! A crystal: the deck as read, and the calculation deck::Build made of it.
extern "C++" struct qcb_solid
{
    qchem::deck::RunSpec                       spec;
    std::unique_ptr<qchem::SolidCalculation>   calc;
    bool                                       attempted = false;   // Result() before any Converge has nothing to say
};

namespace {

using qchem::Calculation;
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

//! A crystal's report run: rundeck brackets Build + Converge with report::Begin/End, and the record its providers write
//! (grids, basis, scf) exists only inside such a bracket.  The bridge splits them across two calls, so each gets its own.
struct RunRecord
{
    explicit RunRecord(const std::string& name) { qchem::report::Begin(name); }
    ~RunRecord() { qchem::report::End(); }
    RunRecord(const RunRecord&) = delete;
    RunRecord& operator=(const RunRecord&) = delete;
};

qchem::SCFParams scf_params(const qcb_request& req)
{
    qchem::SCFParams p;
    p.NMaxIter = size_t(req.maxIter);
    if (req.relaxStart > 0.0) p.StartingRelaxRo = req.relaxStart;
    if (req.relaxMax   > 0.0) p.MaxRelaxRo      = req.relaxMax;
    if (req.relaxMin   > 0.0) p.MinRelaxRo      = req.relaxMin;
    p.UseMOM = req.useMOM != 0;
    if (req.momStart   > 0)   p.MOMStartIter    = req.momStart;
    return p;
}

//! Thrown through the SCF loop when the sink asks to stop (qcb_row_fn's nonzero return); caught at the C boundary.
struct StoppedByUser : std::runtime_error { StoppedByUser() : std::runtime_error("stopped by the user") {} };

//! The facade's observer, forwarding each iteration's row to the C sink: every cell, generically.  One row type serves
//! both facades (Calculation::Observer and SolidCalculation::Observer take the same IterationRow).
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
    }
    return "unknown SCF failure";
}

//! The SCFFailure as one line for Python: why, after how many iterations, the residual, and the facade's details.
std::string describe(const qchem::SCFFailure& f)
{
    char residual[32];
    std::snprintf(residual, sizeof residual, "%.2e", f.residual);
    std::string s = why_name(f.why) + " after " + std::to_string(f.iterations) + " iterations (residual " + residual + ")";
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

void sample(const qcb_calc* h, const Calculation::sf_t& f, int n, double pad, double* out, double* origin, double* spacing)
{
    size_t idx = 0;
    for_grid(h, n, pad, origin, spacing, [&](const rvec3_t& r){ out[idx++] = f(r); });
}

} // anon namespace

extern "C" {

qcb_calc* qcb_make(const qcb_request* req, qcb_row_fn step, qcb_sink* sink, char* err, int errlen)
{
    const HostLocale host;
    try
    {
        auto h = std::make_unique<qcb_calc>();
        h->positions.assign(req->pos3, req->pos3 + 3*req->nat);

        qchem::Molecule mol;
        for (int a = 0; a < req->nat; ++a)
            mol.Insert(new qchem::Atom(req->Z[a], 0.0, rvec3_t(req->pos3[3*a], req->pos3[3*a+1], req->pos3[3*a+2])));

        qchem::CalcOptions opts;
        opts.basis        = req->basis;
        opts.model        = parse_model(req->model);
        opts.pp           = req->pp ? req->pp : "";
        opts.symmetry     = req->symmetry != 0;
        opts.spin         = req->polarized ? qchem::SpinGroup::Polarized : qchem::SpinGroup::UnPolarized;
        opts.multiplicity = req->multiplicity;
        if (req->meshRadial  > 0) opts.mesh.nRadial       = req->meshRadial;
        if (req->meshAngular > 0) opts.mesh.angularDegree = req->meshAngular;
        qchem::AcceleratorOptions acc;
        acc.type = req->accelerator && *req->accelerator ? req->accelerator : "DIIS";
        if (req->diisDepth    > 0)   acc.nProj = req->diisDepth;
        if (req->diisSVTol    > 0.0) acc.svTol = req->diisSVTol;
        if (req->diisEngageFD > 0.0) acc.eMax  = req->diisEngageFD;
        if (req->ladderSwitchAt >= 0.0) acc.ladderSwitchAt = req->ladderSwitchAt;   // 0 = no tail trigger (stall only)

        const auto setup = qchem::MoleculeSetup::Resolve(mol, opts, acc);   // deep-copies mol
        // The observer rides the constructor, so the trace of ITS SCF streams from the first iteration.
        h->params = scf_params(*req);
        h->calc = std::make_unique<Calculation>(setup, h->params,
                                                step ? forward_to(step, sink) : Calculation::Observer{});
        // The sink lives only for this call: detach it before returning (qcb_converge attaches its own).
        h->calc->OnIteration(Calculation::Observer{});
        return h.release();
    }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    return nullptr;
}

void qcb_free(qcb_calc* h) { delete h; }

int qcb_converge(qcb_calc* h, int maxIter, qcb_row_fn step, qcb_sink* sink, char* err, int errlen)
{
    const HostLocale host;
    try
    {
        h->calc->OnIteration(step ? forward_to(step, sink) : Calculation::Observer{});
        qchem::SCFParams p = h->params;
        p.NMaxIter = size_t(maxIter);
        const auto out = h->calc->Converge(p);
        h->calc->OnIteration(Calculation::Observer{});      // the sink lives only for this call
        if (out) return 1;
        write_error(describe(out.Error()).c_str(), err, errlen);
        return 0;
    }
    catch (const StoppedByUser& e) { h->calc->OnIteration(Calculation::Observer{}); write_error(e.what(), err, errlen); return 0; }
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
    for (size_t a = 0; a < st.GetNumAtoms(); ++a) Z_out[a] = int(st[a]->itsZ);
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
        qchem::deck::FromJson(nlohmann::json::parse(deckJson), h->spec);   // strict: an unknown key throws, naming it
        const RunRecord record(h->spec.StructureName() + " build");
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
        const RunRecord record(h->spec.StructureName());
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

int qcb_basis_info(const int* Z, int nat, const double* pos3, const char* basis, double* out4, char* err, int errlen)
{
    const HostLocale host;
    try
    {
        qchem::Molecule mol;
        for (int a = 0; a < nat; ++a)
            mol.Insert(new qchem::Atom(Z[a], 0.0, rvec3_t(pos3[3*a], pos3[3*a+1], pos3[3*a+2])));
        // the facade's defaults (Calculation's BuildBasis with CalcOptions{}): MnD engine, Cartesian, no SALC blocking --
        // a unitary re-blocking would not move S's spectrum
        const std::unique_ptr<qchem::BasisSet::Real_BS> bs(qchem::BasisSet::Gaussian::Factory(
            nlohmann::json{{"basis", basis}, {"engine", "mnd"}, {"angular", "cartesian"}}, &mol));
        double smin = std::numeric_limits<double>::infinity(), smax = 0.0;
        int nblocks = 0;
        for (auto b : bs->Iterate())
        {
            const auto* ov = dynamic_cast<const qchem::BasisSet::Integrals_Overlap<double>*>(b);
            assert(ov && "a molecular orbital block has an overlap");
            rvec_t w; mat_t<double> U;
            qchem::blazem::eigen(ov->Overlap(), w, U);            // ascending
            if (w.size()) { smin = std::min(smin, w[0]); smax = std::max(smax, w[w.size()-1]); }
            ++nblocks;
        }
        out4[0] = double(bs->GetNumFunctions()); out4[1] = nblocks; out4[2] = smin; out4[3] = smax;
        return 1;
    }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    return 0;
}

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
