// File: pybind/qchem_bridge.cpp
//
// Module unit: the global-module-fragment (#includes) de-dups against the
// imported qchem modules' GMFs, so textually including json/STL alongside
// `import qchem.*` no longer triggers std-header ODR errors. The public surface
// is the flat extern "C" API in qcb_api.h (C linkage -> links from nanobind).
module;
#include <nlohmann/json.hpp>
#include <vector>
#include <memory>
#include <algorithm>
#include <utility>
#include <cmath>
#include <cstring>
#include "qcb_api.h"

export module qchem.bridge;

import qchem.Types;                          // rvec3_t (Vector3D<double>), vec3_t
import qchem.Structure;                      // Molecule, Atom, Structure
import qchem.ScalarFunction;                 // ::ScalarFunction<double>
import qchem.BasisSet;                        // Real_BS (GetIrreps)
import qchem.BasisSet.Molecule.Factory;      // BasisSet::Molecule::Factory
import qchem.Hamiltonian;                     // Hamiltonian
import qchem.Hamiltonian.Factory;            // Model, Pol, Factory
import qchem.ElectronConfiguration;          // ElectronConfiguration
import qchem.ElectronConfiguration.Molecule; // Molecule_EC
import qchem.SCFAccelerator;                 // SCFAccelerator
import qchem.SCFAccelerator.Factory;         // SCFAccelerators::Type, Factory
import qchem.SCFIterator;                    // SCFIterator, SCFParams
import qchem.SCFParams;
import qchem.WaveFunction;                   // tWaveFunction (GetChargeDensity/GetOrbitals)
import qchem.Orbitals;                       // Orbitals, Orbital
import qchem.ChargeDensity;                  // DM_CD
import qchem.PeriodicTable;                  // (unused here; symbols mapped in Python)

using qchem::Hamiltonian::Model;
using qchem::Hamiltonian::Pol;
using qchem::SCFIterator::SCFIterator;
using qchem::ChargeDensity::DM_CD;
using qchem::Orbitals::Orbital;
// Public types that moved under qchem:: in the 2026-06 namespace unification.  (Lowercase vocabulary
// like rvec3_t/dcmplx stays global.)  Pulled in by-name rather than `using namespace qchem;` to avoid
// clashing with the SCFIterator/Hamiltonian/Orbitals namespace-vs-class using-declarations above.
using qchem::Molecule;
using qchem::Atom;
using qchem::Structure;
using qchem::ScalarFunction;
using qchem::Real_BS;
using qchem::Spin;
using qchem::ElectronConfiguration;
using qchem::Molecule_EC;
namespace BasisSet = qchem::BasisSet;

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

// Holds the converged SCF; all qchem types stay inside.
struct Calc {
    std::shared_ptr<const Structure> structure;   // also held by the Hamiltonian
    ElectronConfiguration*           ec    = nullptr;  // owned
    BasisSet::Real_BS*               basis = nullptr;  // owned
    SCFIterator*                     scf   = nullptr;   // owns ham + accelerator
    std::vector<int>    numbers;
    std::vector<double> positions;
    ~Calc() { delete scf; delete basis; delete ec; }
};

} // anon namespace

extern "C" {

void* qcb_make(const int* Z, int nat, const double* pos3, const char* basis, int maxiter)
{
    auto* c = new Calc();
    c->numbers.assign(Z, Z + nat);
    c->positions.assign(pos3, pos3 + 3*nat);

    Molecule* mol = new Molecule();
    int Ztot = 0;
    for (int a = 0; a < nat; ++a) {
        mol->Insert(new Atom(Z[a], 0.0, rvec3_t(pos3[3*a], pos3[3*a+1], pos3[3*a+2])));
        Ztot += Z[a];
    }
    std::shared_ptr<Molecule> molp(mol);
    c->structure = molp;
    c->ec    = new Molecule_EC(int(mol->GetNumElectrons()));
    c->basis = BasisSet::Molecule::Factory(nlohmann::json{{"basis", basis}}, mol);

    auto* ham = Factory(Model::HF, Pol::UnPolarized, c->structure);
    nlohmann::json jsacc = {{"NProj", 4}, {"EMax", Ztot*Ztot*0.1/32},
                            {"EMin", 1e-7}, {"SVTol", 5e-9}, {"type", "DIIS"}};
    auto* acc = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);
    c->scf = new SCFIterator(c->basis, c->ec, ham, acc,
                             qchem::ChargeDensity::SeedStrategy::Default, mol);
    //          NMaxIter MinΔρ MinΔE MinVirial MinFD relax mergeTol verbose
    c->scf->Iterate({size_t(maxiter), 1e-4, 1e-7, 1e-13, 1e-5, 1.0, 1e-4, false});
    return c;
}

void   qcb_free(void* h)        { delete static_cast<Calc*>(h); }
double qcb_energy(void* h)      { return static_cast<Calc*>(h)->scf->GetEnergy().GetTotalEnergy(); }
int    qcb_natoms(void* h)      { return int(static_cast<Calc*>(h)->numbers.size()); }

void qcb_atoms(void* h, int* Z_out, double* xyz_out)
{
    auto* c = static_cast<Calc*>(h);
    std::memcpy(Z_out, c->numbers.data(), c->numbers.size()*sizeof(int));
    std::memcpy(xyz_out, c->positions.data(), c->positions.size()*sizeof(double));
}

void qcb_density(void* h, int n, double pad, double* out, double* origin, double* spacing)
{
    auto* c = static_cast<Calc*>(h);
    rvec3_t lo, hi; bbox(c->positions, pad, lo, hi);
    std::unique_ptr<DM_CD> cd(c->scf->GetWaveFunction()->GetChargeDensity());  // caller owns
    sample_scalar(*cd, lo, hi, n, out, origin, spacing);
}

void qcb_gradient(void* h, int n, double pad, double* out, double* origin, double* spacing)
{
    auto* c = static_cast<Calc*>(h);
    rvec3_t lo, hi; bbox(c->positions, pad, lo, hi);
    std::unique_ptr<DM_CD> cd(c->scf->GetWaveFunction()->GetChargeDensity());
    rvec3_t sp((hi.x-lo.x)/(n-1), (hi.y-lo.y)/(n-1), (hi.z-lo.z)/(n-1));
    origin[0]=lo.x; origin[1]=lo.y; origin[2]=lo.z;
    spacing[0]=sp.x; spacing[1]=sp.y; spacing[2]=sp.z;
    size_t idx = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                vec3_t<double> g = cd->Gradient(rvec3_t(lo.x+sp.x*i, lo.y+sp.y*j, lo.z+sp.z*k));
                out[idx++] = g.x; out[idx++] = g.y; out[idx++] = g.z;
            }
}

void qcb_orbital(void* h, int index, int n, double pad,
                 double* out, double* origin, double* spacing, int* is_signed)
{
    auto* c = static_cast<Calc*>(h);
    const auto* wf = c->scf->GetWaveFunction();
    std::vector<std::pair<double, const ScalarFunction<double>*>> occ;
    for (const auto& irr : c->basis->GetIrreps(Spin::None))
        for (const Orbital* o : wf->GetOrbitals(irr)->Iterate())
            if (o->IsOccupied()) {
                const auto* phi = dynamic_cast<const ScalarFunction<double>*>(o);
                if (phi) occ.push_back({o->GetEigenEnergy(), phi});
            }
    std::sort(occ.begin(), occ.end(), [](const auto& a, const auto& b){ return a.first > b.first; });
    *is_signed = 1;
    rvec3_t lo, hi; bbox(c->positions, pad, lo, hi);
    if (occ.empty() || index < 0 || size_t(index) >= occ.size()) {
        std::memset(out, 0, size_t(n)*n*n*sizeof(double));
        origin[0]=lo.x; origin[1]=lo.y; origin[2]=lo.z;
        spacing[0]=spacing[1]=spacing[2]=0;
        return;
    }
    sample_scalar(*occ[index].second, lo, hi, n, out, origin, spacing);
}

} // extern "C"
