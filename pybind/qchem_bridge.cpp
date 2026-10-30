// File: pybind/qchem_bridge.cpp
//
// Module unit: the global-module-fragment (#includes) de-dups against the imported qchem modules' GMFs, so
// textually including STL alongside `import qchem.*` does not trigger std-header ODR errors.  The public
// surface is the flat extern "C" API in qcb_api.h (C linkage -> links from nanobind).
//
// The bridge holds NO recipe of its own: it fills CalcOptions / AcceleratorOptions / SCFParams, calls
// MoleculeSetup::Resolve, and owns one qchem::Calculation.  Everything it answers comes from the facade.
module;
#include <algorithm>
#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <vector>
#include "qcb_api.h"

export module qchem.bridge;

import qchem.Types;          // rvec3_t, vec3_t
import qchem.Structure;      // Molecule, Atom
import qchem.SCFParams;      // SCFParams
import qchem.SCFIterator;    // SCFProgress
import qchem.Calculation;    // Calculation, MoleculeSetup, CalcOptions, AcceleratorOptions

// The opaque handle: the facade plus the request's geometry (for the sampling box).  extern "C++" attaches it to the
// GLOBAL module, where qcb_api.h (in the GMF) forward-declares it -- a purview definition would be a different entity.
extern "C++" struct qcb_calc
{
    std::unique_ptr<qchem::Calculation> calc;
    std::vector<double>                 positions;   // flat 3N bohr, as asked
};

namespace {

using qchem::Calculation;

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

qchem::SCFParams scf_params(int maxIter)
{
    qchem::SCFParams p;
    p.NMaxIter = size_t(maxIter);
    return p;
}

//! The facade's observer, forwarding each iteration to the C sink.
Calculation::Observer forward_to(qcb_step_fn step, qcb_sink* sink)
{
    return [step, sink](const qchem::SCFIterator::SCFProgress& p)
    {
        const qcb_step s{int(p.iteration), p.energy, p.dE, p.commutator, p.drho};
        step(sink, &s);
    };
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

qcb_calc* qcb_make(const qcb_request* req, qcb_step_fn step, qcb_sink* sink, char* err, int errlen)
{
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
        opts.multiplicity = req->multiplicity;
        qchem::AcceleratorOptions acc;
        acc.type = req->accelerator && *req->accelerator ? req->accelerator : "DIIS";

        const auto setup = qchem::MoleculeSetup::Resolve(mol, opts, acc);   // deep-copies mol
        // The observer rides the constructor, so the trace of ITS SCF streams from the first iteration.
        h->calc = std::make_unique<Calculation>(setup, scf_params(req->maxIter),
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

int qcb_converge(qcb_calc* h, int maxIter, qcb_step_fn step, qcb_sink* sink, char* err, int errlen)
{
    try
    {
        h->calc->OnIteration(step ? forward_to(step, sink) : Calculation::Observer{});
        return h->calc->Converge(scf_params(maxIter)) ? 1 : 0;
    }
    catch (const std::exception& e) { write_error(e.what(), err, errlen); }
    catch (...)                     { write_error("qchem_py: unknown failure", err, errlen); }
    return -1;
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

} // extern "C"
