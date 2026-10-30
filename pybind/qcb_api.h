// File: pybind/qcb_api.h
//
// Flat C ABI between the two worlds:
//   qchem_bridge.cpp -- a C++20 *module unit* that `import`s qchem (its #includes live in a
//                       global-module-fragment, so they de-dup against the imported modules' GMFs --
//                       no <ranges> ODR clash).  It runs MoleculeSetup::Resolve + qchem::Calculation.
//   qchem_py.cpp     -- pure nanobind; includes ONLY this header, never a module.
//
// HANDLES ARE OPAQUE STRUCT POINTERS, never void*.  `qcb_calc` is defined only inside the bridge (it owns the
// Calculation); `qcb_sink` is defined only inside the nanobind TU (it holds the Python callable the SCF trace
// streams to).  Each side sees the other's type as incomplete, so the compiler still checks every pointer.
//
// extern "C" gives every function C language linkage, so the symbols are plain (unmangled, not attached to the
// module) and link cleanly from the nanobind TU.  Buffers are caller-allocated and fixed-size (the caller knows
// n), so no STL or ownership crosses.  Nothing here throws: a failure returns null / -1 and writes a message
// into the caller's `err` buffer.
#pragma once
#ifdef __cplusplus
extern "C" {
#endif

typedef struct qcb_calc qcb_calc;   //!< a converged molecular calculation (bridge-owned)
typedef struct qcb_sink qcb_sink;   //!< the caller's per-iteration sink (nanobind-owned)

//! One SCF iteration, as the facade's observer reports it.
typedef struct qcb_step
{
    int    iteration;
    double energy;      //!< total energy (hartree)
    double dE;          //!< |E_n - E_{n-1}|
    double commutator;  //!< [F,D]
    double drho;        //!< density change
} qcb_step;
typedef void (*qcb_step_fn)(qcb_sink*, const qcb_step*);

//! What the caller asks for: a molecule (Z + flat 3N bohr positions) and the CalcOptions the GUI exposes.
typedef struct qcb_request
{
    const int*    Z;
    int           nat;
    const double* pos3;
    const char*   basis;          //!< "sto-3g" | "dzvp" | ...
    const char*   model;          //!< "HF" | "LDA" | "Xalpha"
    const char*   pp;             //!< pseudopotential set name, "" = all-electron
    const char*   accelerator;    //!< "DIIS" | "GDM" | "Ladder"
    int           symmetry;       //!< point-group SALC blocking
    int           multiplicity;   //!< 2S+1, 0 = minimal spin
    int           maxIter;
} qcb_request;

//! Resolve + converge.  `step`/`sink` may be null; when given, the trace streams live.  Null on failure (`err`).
qcb_calc* qcb_make(const qcb_request* req, qcb_step_fn step, qcb_sink* sink, char* err, int errlen);
void      qcb_free(qcb_calc* h);
//! Re-run the SCF from the seed, streaming each iteration to `step`.  Returns converged (1) / not (0) / failed (-1).
int       qcb_converge(qcb_calc* h, int maxIter, qcb_step_fn step, qcb_sink* sink, char* err, int errlen);

double qcb_energy    (const qcb_calc* h);
int    qcb_converged (const qcb_calc* h);
int    qcb_iterations(const qcb_calc* h);
int    qcb_natoms    (const qcb_calc* h);
void   qcb_atoms     (const qcb_calc* h, int* Z_out, double* xyz_out);   // Z_out[nat], xyz_out[3*nat]

//! Every MO, virtuals included, ascending in energy: n = qcb_num_levels of them.  spin: -1 down, 0 none, +1 up.
//! `label` is nlabel chars per level.
int    qcb_num_levels(const qcb_calc* h);
void   qcb_levels(const qcb_calc* h, double* e, double* occ, int* spin, char* label, int nlabel);

// Sample fields on an n^3 grid padded `pad` bohr past the atoms. The caller preallocates `out` (n*n*n for
// scalars, n*n*n*3 for the vector field, C-order) plus origin[3] and spacing[3].  `index` is a LEVEL index
// (qcb_levels' order, 0 = lowest).  qcb_orbital returns 0 for an index out of range.  The spin density of an
// unpolarized run is the zero field.
void qcb_density     (const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing);
void qcb_spin_density(const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing);
int  qcb_orbital     (const qcb_calc* h, int index, int n, double pad, double* out, double* origin, double* spacing);
void qcb_gradient(const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing);

#ifdef __cplusplus
}
#endif
