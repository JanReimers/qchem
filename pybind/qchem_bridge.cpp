// File: pybind/qchem_bridge.cpp
//
// Module unit bridging qchem -> the flat extern "C" API (qcb_api.h). Now built on
// the qchem::Calculation front door: one `import qchem;` (the umbrella) gives the
// facade + Molecule/Atom + ScalarFunction sampling, so this file is just "build a
// Calculation, sample its Density()/Orbital()/Gradient() onto grids". The old
// hand-assembly (EC + basis Factory + Hamiltonian + accelerator + SCFIterator)
// now lives once, in the lib. extern "C" gives plain C symbols -> links from the
// nanobind TU (qchem_py.cpp), which never touches a C++20 module.
module;
#include <vector>
#include <memory>
#include <cstring>
#include "qcb_api.h"

export module qchem.bridge;
import qchem;   // umbrella: Calculation, CalcOptions, Molecule, Atom, ScalarFunction, rvec3_t, SCFProgress

using qchem::Calculation;
using qchem::CalcOptions;
using qchem::Molecule;
using qchem::Atom;
using qchem::ScalarFunction;

namespace {

void bbox(const std::vector<double>& pos, double pad, rvec3_t& lo, rvec3_t& hi)
{
    double mn[3] = { 1e30,  1e30,  1e30}, mx[3] = {-1e30, -1e30, -1e30};
    for (size_t a = 0; a < pos.size()/3; ++a)
        for (int d = 0; d < 3; ++d) {
            mn[d] = std::min(mn[d], pos[3*a+d]);
            mx[d] = std::max(mx[d], pos[3*a+d]);
        }
    lo = rvec3_t(mn[0]-pad, mn[1]-pad, mn[2]-pad);
    hi = rvec3_t(mx[0]+pad, mx[1]+pad, mx[2]+pad);
}

// Fill out[n^3] (C-order, k fastest) with f(r); write origin/spacing.
void sample_scalar(const ScalarFunction<double>& f, rvec3_t lo, rvec3_t hi, int n,
                   double* out, double* origin, double* spacing)
{
    rvec3_t sp((hi.x-lo.x)/(n-1), (hi.y-lo.y)/(n-1), (hi.z-lo.z)/(n-1));
    origin[0]=lo.x; origin[1]=lo.y; origin[2]=lo.z;
    spacing[0]=sp.x; spacing[1]=sp.y; spacing[2]=sp.z;
    size_t idx = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                out[idx++] = f(rvec3_t(lo.x+sp.x*i, lo.y+sp.y*j, lo.z+sp.z*k));
}

// The converged calculation + the geometry we were handed (for bbox / atom readout).
struct Calc {
    std::unique_ptr<Calculation> calc;
    std::vector<int>    numbers;
    std::vector<double> positions;
    int                 maxiter = 20;
};

} // anon namespace

extern "C" {

void* qcb_make(const int* Z, int nat, const double* pos3, const char* basis, int maxiter)
{
    auto* c = new Calc();
    c->numbers.assign(Z, Z + nat);
    c->positions.assign(pos3, pos3 + 3*nat);
    c->maxiter = maxiter;

    Molecule mol;   // local; Calculation deep-copies it
    for (int a = 0; a < nat; ++a)
        mol.Insert(new Atom(Z[a], 0.0, rvec3_t(pos3[3*a], pos3[3*a+1], pos3[3*a+2])));

    c->calc = std::make_unique<Calculation>(mol, CalcOptions{.basis = basis});  // build + converge
    if (maxiter != 20)
        c->calc->Converge({.NMaxIter = size_t(maxiter)});                        // honor a custom cap
    return c;
}

int qcb_run_scf(void* h, qcb_scf_cb cb, void* user)
{
    auto* c = static_cast<Calc*>(h);
    c->calc->OnIteration([cb, user](const qchem::SCFIterator::SCFProgress& p) {
        cb(user, int(p.iteration), p.energy, p.dE, p.commutator, p.drho);
    });
    c->calc->Converge({.NMaxIter = size_t(c->maxiter)});   // re-run from the seed, streaming
    return int(c->calc->IterationCount());
}

void   qcb_free(void* h)   { delete static_cast<Calc*>(h); }
double qcb_energy(void* h) { return static_cast<Calc*>(h)->calc->Energy(); }
int    qcb_natoms(void* h) { return int(static_cast<Calc*>(h)->numbers.size()); }

void qcb_atoms(void* h, int* Z_out, double* xyz_out)
{
    auto* c = static_cast<Calc*>(h);
    std::memcpy(Z_out,   c->numbers.data(),   c->numbers.size()  * sizeof(int));
    std::memcpy(xyz_out, c->positions.data(), c->positions.size()* sizeof(double));
}

void qcb_density(void* h, int n, double pad, double* out, double* origin, double* spacing)
{
    auto* c = static_cast<Calc*>(h);
    rvec3_t lo, hi; bbox(c->positions, pad, lo, hi);
    sample_scalar(c->calc->Density(), lo, hi, n, out, origin, spacing);
}

void qcb_gradient(void* h, int n, double pad, double* out, double* origin, double* spacing)
{
    auto* c = static_cast<Calc*>(h);
    rvec3_t lo, hi; bbox(c->positions, pad, lo, hi);
    const ScalarFunction<double>& rho = c->calc->Density();
    rvec3_t sp((hi.x-lo.x)/(n-1), (hi.y-lo.y)/(n-1), (hi.z-lo.z)/(n-1));
    origin[0]=lo.x; origin[1]=lo.y; origin[2]=lo.z;
    spacing[0]=sp.x; spacing[1]=sp.y; spacing[2]=sp.z;
    size_t idx = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                vec3_t<double> g = rho.Gradient(rvec3_t(lo.x+sp.x*i, lo.y+sp.y*j, lo.z+sp.z*k));
                out[idx++] = g.x; out[idx++] = g.y; out[idx++] = g.z;
            }
}

void qcb_orbital(void* h, int index, int n, double pad,
                 double* out, double* origin, double* spacing, int* is_signed)
{
    auto* c = static_cast<Calc*>(h);
    *is_signed = 1;
    rvec3_t lo, hi; bbox(c->positions, pad, lo, hi);
    const size_t nocc = c->calc->NumOccupied();
    if (index < 0 || size_t(index) >= nocc) {        // out of range -> empty field
        std::memset(out, 0, size_t(n)*n*n*sizeof(double));
        origin[0]=lo.x; origin[1]=lo.y; origin[2]=lo.z;
        spacing[0]=spacing[1]=spacing[2]=0;
        return;
    }
    // index 0 = HOMO; Calculation::Orbital(i) is energy-ascending (0 = lowest).
    sample_scalar(c->calc->Orbital(nocc - 1 - size_t(index)), lo, hi, n, out, origin, spacing);
}

} // extern "C"
