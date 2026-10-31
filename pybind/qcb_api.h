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
typedef struct qcb_solid qcb_solid; //!< a periodic calculation built from a deck (bridge-owned)

//! One cell of an SCF iteration row: one column of the console trace (qchem::SCFIterator::IterationCell), GENERIC -- a
//! column the lib adds to a layout reaches Python with no ABI change.  The strings are UTF-8 ("Δρ") and, like the
//! cells array, valid ONLY for the duration of the callback.
typedef struct qcb_cell
{
    const char* title;      //!< the console header ("Etotal", "[F,D]", "Δρ", "m_site", ...)
    double      value;      //!< the number a plot draws; NaN when the cell has none (a "----", a text-only column)
    const char* text;       //!< the cell exactly as the console prints it, padding stripped
    double      tolerance;  //!< the gate under the header; NaN when the column gates nothing
    int         logScale;   //!< 1: a convergence measure (plot |value| on a log axis); 0: linear
} qcb_cell;

//! One SCF iteration: its cells in display order ("#" excluded) and the total energy (the row's EnergyBreakdown total).
typedef struct qcb_row
{
    int             iteration;
    double          energy;   //!< total energy (hartree)
    int             ncells;
    const qcb_cell* cells;
} qcb_row;
typedef void (*qcb_row_fn)(qcb_sink*, const qcb_row*);

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
    int           polarized;      //!< spin-native two-channel run (SpinGroup::Polarized) even for a closed shell
    int           multiplicity;   //!< 2S+1, 0 = minimal spin; > 1 promotes to polarized anyway
    int           maxIter;
    int           meshRadial;     //!< XC mesh radial points per atom, 0 = the facade's default (CalcOptions::mesh)
    int           meshAngular;    //!< XC mesh Lebedev degree, 0 = the facade's default
} qcb_request;

//! Resolve + converge.  `step`/`sink` may be null; when given, the trace streams live.  Null on failure (`err`).
qcb_calc* qcb_make(const qcb_request* req, qcb_row_fn step, qcb_sink* sink, char* err, int errlen);
void      qcb_free(qcb_calc* h);
//! Re-run the SCF from the seed, streaming each iteration to `step`.  Returns converged (1) / not (0) / failed (-1).
//! On 0 the SCFFailure's reason is written to `err` (as qcb_failure does).
int       qcb_converge(qcb_calc* h, int maxIter, qcb_row_fn step, qcb_sink* sink, char* err, int errlen);
//! Why the last SCF did not converge (SCFFailure: why, iterations, residual, details) written to `out`; returns 0 and
//! leaves `out` empty when it converged.
int       qcb_failure(const qcb_calc* h, char* out, int outlen);

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

// -- periodic runs: the INPUT DECK is the seam (qchem::deck::Build).  The GUI translates its own run spec into deck JSON
//    (the same file `rundeck` runs), so the bridge restates none of the deck's resolve logic and every GUI crystal run is a
//    deck `rundeck` reproduces.

//! Parse `deckJson` (schema 1, src/Calculation/Deck.C) and BUILD its calculation, stopping before the SCF (deck::Build:
//! resolve, lattice, basis, SolidSetup).  Null on a refusal -- a bad deck, an unknown key, a setup the library refuses --
//! with the library's own text in `err`.
qcb_solid* qcb_solid_make(const char* deckJson, char* err, int errlen);
void       qcb_solid_free(qcb_solid* h);
//! Converge the deck's schedule (deck::ScheduleToRun), streaming each iteration to `step`.  Returns converged (1) / not
//! (0, the SCFFailure's reason in `err`) / failed (-1, the exception's text in `err`).
int        qcb_solid_converge(qcb_solid* h, qcb_row_fn step, qcb_sink* sink, char* err, int errlen);
//! As qcb_failure, for the last qcb_solid_converge (before any: "not converged yet").
int        qcb_solid_failure(const qcb_solid* h, char* out, int outlen);
//! The converged total energy (hartree); NaN when the last SCF did not converge -- a failed run's last iterate is not an answer.
double     qcb_solid_energy    (const qcb_solid* h);
int        qcb_solid_iterations(const qcb_solid* h);
//! The resolved deck (every default filled in, Nelec / species derived): what `rundeck` would record.  Returns the length
//! it needs (excluding the NUL); writes at most outlen-1 chars.
int        qcb_solid_resolved(const qcb_solid* h, char* out, int outlen);

#ifdef __cplusplus
}
#endif
