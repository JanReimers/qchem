// File: pybind/qcb_api.h
//
// Flat C ABI between the two worlds:
//   qchem_bridge.cpp -- a C++20 *module unit* that `import`s qchem (its #includes
//                       live in a global-module-fragment, so they de-dup against
//                       the imported modules' GMFs -- no <ranges> ODR clash).
//   qchem_py.cpp     -- pure nanobind; includes ONLY this header, never a module.
//
// extern "C" gives every function C language linkage, so the symbols are plain
// (unmangled, not attached to the module) and link cleanly from the nanobind TU.
// All buffers are caller-allocated and fixed-size (the caller knows n), so no
// STL or ownership crosses the boundary.
#pragma once
#ifdef __cplusplus
extern "C" {
#endif

// Build a molecule (atomic numbers + flat 3N bohr positions) and converge an SCF.
// method: "HF" (default) | "LDA" | "Xalpha" -> CalcOptions.model.
void*  qcb_make(const int* Z, int nat, const double* pos3,
                const char* basis, const char* method, int maxiter);
void   qcb_free(void* h);

double qcb_energy(void* h);
int    qcb_natoms(void* h);
void   qcb_atoms (void* h, int* Z_out, double* xyz_out);   // Z_out[nat], xyz_out[3*nat]

// Re-run the SCF from the seed, calling `cb` once per iteration (live convergence).
// Leaves the calculator in a freshly-converged state. Returns the iteration count.
typedef void (*qcb_scf_cb)(void* user, int iter, double E, double dE,
                           double commutator, double drho);
int qcb_run_scf(void* h, qcb_scf_cb cb, void* user);

// Sample fields on an n^3 grid padded `pad` bohr past the atoms. The caller
// preallocates `out` (n*n*n for scalars, n*n*n*3 for the vector field, C-order)
// plus origin[3] and spacing[3]. Density/orbital/gradient all sample the same
// converged state.
void qcb_density (void* h, int n, double pad, double* out, double* origin, double* spacing);
void qcb_orbital (void* h, int index, int n, double pad,
                  double* out, double* origin, double* spacing, int* is_signed);
void qcb_gradient(void* h, int n, double pad, double* out, double* origin, double* spacing);

#ifdef __cplusplus
}
#endif
