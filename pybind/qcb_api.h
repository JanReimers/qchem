// File: pybind/qcb_api.h
//
// Flat C ABI between the two worlds:
//   qchem_bridge.cpp -- a C++20 *module unit* that `import`s qchem (its #includes live in a
//                       global-module-fragment, so they de-dup against the imported modules' GMFs --
//                       no <ranges> ODR clash).  It runs MoleculeSetup::Resolve + qchem::MoleculeCalculation.
//   qchem_py.cpp     -- pure nanobind; includes ONLY this header, never a module.
//
// HANDLES ARE OPAQUE STRUCT POINTERS, never void*.  `qcb_calc` is defined only inside the bridge (it owns the
// MoleculeCalculation); `qcb_sink` is defined only inside the nanobind TU (it holds the Python callable the SCF trace
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
//! Returns 0 to go on, nonzero to STOP the SCF (the user's Stop button): the bridge then abandons the solve and reports
//! "stopped by the user" (qcb_make: null; the converge calls: 0).
typedef int (*qcb_row_fn)(qcb_sink*, const qcb_row*);

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
    double        relaxStart;     //!< linear D-mixing's starting step (SCFParams::StartingRelaxRo), 0 = the lib's default
    double        relaxMax;       //!< its adaptive ceiling (SCFParams::MaxRelaxRo), 0 = the lib's default
    int           diisDepth;      //!< DIIS history (AcceleratorOptions::nProj), 0 = the lib's default
    double        diisSVTol;      //!< DIIS singular-value floor (AcceleratorOptions::svTol), 0 = the lib's default
    double        diisEngageFD;   //!< DIIS engages below this [F,D] (AcceleratorOptions::eMax), 0 = the lib's Z-scaled default
    int           useMOM;         //!< Maximum Overlap Method (SCFParams::UseMOM)
    int           momStart;       //!< its delayed reference capture (SCFParams::MOMStartIter), 0 = the lib's default
    double        relaxMin;       //!< the adaptive step's floor (SCFParams::MinRelaxRo), 0 = the lib's default
    double        ladderSwitchAt; //!< the Ladder's DIIS -> GDM hand-off at |ΔE/E| below this (AcceleratorOptions::ladderSwitchAt),
                                  //!< 0 = stall only, < 0 = the lib's default (1e-6)
    //! The SCF's convergence gates (SCFParams' std::optional Gates): > 0 the gate, 0 = OFF, < 0 = the lib's default.
    double        gateFD;         //!< [F,D], the worst irrep's orbital gradient (MinFD; default 1e-5)
    double        gateDrho;       //!< Δρ per electron (MinΔρPerElectron; default 1e-4)
    double        gateDE;         //!< |ΔE/E| (MinΔE; default off)
    double        gateDFD;        //!< |Δ[F,D]| (MinΔFD; default off)
    const char*   mixer;          //!< SCFParams::Mixing by name: "Auto" | "Linear" | "Kerker" | "Pulay" ("" = Auto); the facade refuses
                                  //!< what a molecule cannot mix (Kerker needs G space)
    double        smearingkT;     //!< Fermi-Dirac smearing kT, Ha (SCFParams::SmearingkT), 0 = off
    const char*   basisFile;      //!< a .bsd path used instead of the shipped file `basis` names (CalcOptions::basisFile); null/"" = none
} qcb_request;

//! Resolve + converge.  `step`/`sink` may be null; when given, the trace streams live.  Null on failure (`err`).
qcb_calc* qcb_make(const qcb_request* req, qcb_row_fn step, qcb_sink* sink, char* err, int errlen);
//! qcb_make's general form.  `path` (null or "" = a fresh start): the SCF starts from that saved state
//! (MoleculeCalculation::Restart; "qchem-state 2.0"): an EXACT RESUME when nothing differs, else a WARM START naming each
//! difference (qcb_restarted_from).  `buildOnly` != 0: BUILD the whole graph and stop before iteration 1
//! (MoleculeCalculation::Prepare / PrepareRestart), filed as the report run "<ID> build" (basis.perIrrep with the overlap's
//! conditioning, every mesh); the first qcb_converge then runs the SCF on that graph.  Null on failure; a REFUSAL of the state
//! (its kind, basis, electrons or blocks differ, or the file is unreadable) also sets *refused = 1, so the caller can tell it
//! from a failure and decide on a fresh run.
qcb_calc* qcb_start(const qcb_request* req, const char* path, int buildOnly, qcb_row_fn step, qcb_sink* sink, int* refused,
                    char* err, int errlen);
void      qcb_free(qcb_calc* h);
//! SAVE the run's state to `path` (HDF5 "qchem-state 2.0", converged or not: a restart says which), filed as its own report run
//! "<ID> state save".  1 on success, 0 with `err`.
int       qcb_save_state(const qcb_calc* h, const char* path, char* err, int errlen);
//! WHAT A RESTART OF THIS RUN FROM `path` WOULD DO, without doing it (Calculation::JudgeState), as JSON text:
//! {"exact": B, "warm": [line, ...]} or {"refused": {"why": unreadable|format|kind|mismatch|realness, "details": S}}.
//! Returns the length; call with outlen 0 to size it.
int       qcb_judge_state(const qcb_calc* h, const char* path, char* out, int outlen);
//! What this run's restart did (Calculation::RestartedFrom), as JSON text: {path, exact, warm: [...], maxDisplacement (bohr),
//! electrons (the restored density's Tr DS on this basis)}, or "null" for a fresh run.  Returns the length.
int       qcb_restarted_from(const qcb_calc* h, char* out, int outlen);
//! Run the SCF (after a buildOnly qcb_start: on its graph; else again from the seed), streaming each iteration to `step`.  Returns converged (1) / not (0) / failed (-1).
//! On 0 the SCFFailure's reason is written to `err` (as qcb_failure does).
int       qcb_converge(qcb_calc* h, int maxIter, qcb_row_fn step, qcb_sink* sink, char* err, int errlen);
//! Why the last SCF did not converge (SCFFailure: why, iterations, details) written to `out`; returns 0 and
//! leaves `out` empty when it converged.
int       qcb_failure(const qcb_calc* h, char* out, int outlen);

double qcb_energy    (const qcb_calc* h);
int    qcb_converged (const qcb_calc* h);
int    qcb_iterations(const qcb_calc* h);
int    qcb_natoms    (const qcb_calc* h);
void   qcb_atoms     (const qcb_calc* h, int* Z_out, double* xyz_out);   // Z_out[nat], xyz_out[3*nat]

//! EVERY LEVEL, virtuals included, ascending in energy (Calculation::MakeLevels, the structure-neutral face; qcb_solid_levels is
//! the crystal's door to the SAME call): JSON [{"e","occ","spin" (-1/0/+1),"k" (fractional),"weight","label"}], or
//! {"error": why} before the first converge.  Returns the JSON's length (call with outlen 0 to size it).
int    qcb_levels(const qcb_calc* h, char* out, int outlen);

// Sample fields on an n^3 grid padded `pad` bohr past the atoms. The caller preallocates `out` (n*n*n for
// scalars, n*n*n*3 for the vector field, C-order) plus origin[3] and spacing[3].  `index` is a LEVEL index
// (qcb_levels' order, 0 = lowest).  qcb_orbital returns 0 for an index out of range.  The spin density of an
// unpolarized run is the zero field.
void qcb_density     (const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing);
void qcb_spin_density(const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing);
int  qcb_orbital     (const qcb_calc* h, int index, int n, double pad, double* out, double* origin, double* spacing);
void qcb_gradient(const qcb_calc* h, int n, double pad, double* out, double* origin, double* spacing);
//! The CURRENT ITERATE's density (spin != 0: its spin density, the zero field when unpolarized), sampled as qcb_density does:
//! the wave function's latest diagonalisation, built fresh per call (MoleculeCalculation::MakeIterateDensity).  Call it ONLY
//! from the row callback (`step`, on the compute thread between iterations) or between converges -- never from another thread
//! while an SCF runs.  1 on success, 0 with `err` (no iteration yet).
int  qcb_iterate_density(const qcb_calc* h, int spin, int n, double pad, double* out, double* origin, double* spacing,
                         char* err, int errlen);

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

//! A crystal's density (spin != 0: its magnetization m = ρ↑−ρ↓; the zero field when unpolarized) on the Cartesian box the
//! caller states: n3[0] x n3[1] x n3[2] points (C-order, the last fastest) at origin + (i, j, k) * spacing, bohr -- the unit cell's
//! bounding box or a supercell's, the caller's choice (ρ is periodic, so any box is valid).  iterate == 0: the CONVERGED density
//! (Converged::Density / SpinDensity; refused, saying so, when the last SCF did not converge); iterate != 0: the CURRENT iterate's
//! (SolidCalculation::MakeIterateDensity -- only from the row callback, or between converges).  `out` holds the product of n3.
//! 1 on success, 0 with `err`.
int        qcb_solid_field(const qcb_solid* h, int spin, int iterate, const int* n3, const double* origin, const double* spacing,
                           double* out, char* err, int errlen);
//! As qcb_judge_state / qcb_restarted_from, for a crystal (the deck saves and restarts through state.save / state.restartFrom).
int        qcb_solid_judge_state(const qcb_solid* h, const char* path, char* out, int outlen);
int        qcb_solid_restarted_from(const qcb_solid* h, char* out, int outlen);
int        qcb_solid_levels(const qcb_solid* h, char* out, int outlen);   //!< as qcb_levels

// -- symmetry DETECTED from a structure (the lib's own detectors; nothing is imposed here) -------------------------------

//! The library's run record (report::GlobalReport) as JSON text: every run since the last clear, keyed "<name>@<start>", each
//! holding the sections its providers emitted (basis, grids, scf, ...).  Returns the length; call with outlen 0 to size it.
int  qcb_report(char* out, int outlen);
//! The OPEN run's record so far (report::CurrentRunReport) as JSON text, "" when no run is open.  For a row callback: it
//! runs on the thread that writes the report, so the snapshot is race-free (what setup filed is there by iteration 1).
int  qcb_current_report(char* out, int outlen);
//! Drop the accumulated run record (a no-op while a run is open).
void qcb_clear_report(void);

//! The molecular point group of `nat` points (`species` labels, flat 3N bohr `pos3`) within `tol` bohr: the full Schoenflies
//! symbol into `symbol`, the abelian subgroup the SCF blocks by into `abelian`.  Returns the group order (0 for C∞v / D∞h),
//! -1 on failure (`symbol` holds the reason).
int qcb_point_group(const int* species, int nat, const double* pos3, double tol,
                    char* symbol, int symlen, char* abelian, int ablen);
//! The space group of a cell (`lattice9`: the three lattice vectors as ROWS, bohr) and a fractional basis (`species`, flat 3N
//! `frac3`) within `tol`: writes up to `maxops` linear parts W (lattice coordinates, row-major 3x3 each) into `W9` and whether
//! every tau is 0 into `symmorphic`.  Returns the group order (-1 on failure, the reason in `err`).
int qcb_space_group(const double* lattice9, const int* species, int nat, const double* frac3, double tol,
                    double* W9, int maxops, int* symmorphic, char* err, int errlen);

// -- crystal FILES and the MATERIALS layers (qchem.CrystalFiles, StructureData's layers; OpenWork §2, ruled 2026-11-09) --------

//! One materials request as JSON text, its answer as JSON text (`{"error": "..."}` on a refusal, in the library's words).  The
//! library reads no environment: the caller passes every library it layers over the shipped one, as
//! `"libraries": [{"layer": "user"|"workspace", "path": P}, ...]` (the GUI resolves ~/.solid2/materials.json itself):
//!   {"op":"list","libraries":L}                          -> {"materials":[{name, layer: shipped|user|workspace, file, shadowed, entry, lattice_bohr}]}
//!   {"op":"import","path":P}                             -> {entry, name, warnings, lattice_bohr}   (CIF / POSCAR / extended xyz)
//!   {"op":"export","entry":E,"format":F,"name":N}        -> {text}   (F: cif | poscar | xyz)
//!   {"op":"save","into":{layer,path},"name":N,"entry":E,"replace":B}  -> {path}   (a later layer SHADOWS an earlier one: a list's
//!                                                          `shadowed` entry is one a later layer's same name hides)
//!   {"op":"cell","entry":E}                              -> {lattice_bohr}   (rows a, b, c; a Bravais entry's vectors)
//! TWO CALLS: with `out` null the request RUNS and its answer is kept (the return is its length); then a call with a buffer
//! copies that answer (`req` ignored).  So a save runs once.
int qcb_materials(const char* req, char* out, int outlen);

//! A Gaussian94 basis from outside (the Basis Set Exchange's "gaussian94" text) TRANSLATED into the BasisSetData layout
//! (BasisSet::Gaussian::ImportGaussian94): {"text": T, "title": S, "provenance": [lines]} -> {bsd, elements, maxL, notes}, or
//! {"error": why} when the library cannot represent it (an ECP block, a scaled shell, ...).  The text is parsed as data only.
//! The same two calls as qcb_materials.
int qcb_basis_import(const char* req, char* out, int outlen);

//! A named pseudopotential set's choice (Pseudopotential::PseudopotentialSet, the ONE reader of pseudopotentials.json and its rule):
//! {"set": S, "elements": [El, ...]} -> {set, functional, valence: {El: q}, fitted: {El: label}} -- q per element, the set's exception if it lists one, else
//! the smallest q the GTH database holds at the set's functional; `fitted` the functional each resolved potential was FITTED with
//! (PseudopotentialSet::FittedFunctional, a display label such as "PADE" -- compare functionals through the library, never by this name);
//! or {"error"} (an unknown set, an element the database lacks).
//! The same two calls as qcb_materials.
int qcb_pseudopotentials(const char* req, char* out, int outlen);

//! A RUN's thread counts as a scope (qchem::ScopedWorkerThreads; OpenWork §2 "Per-run thread count", user 2026-11-12): while open,
//! the OpenMP regions run `openmp` threads and the BLAS pool `blas` (each when > 0; 0 = the process's own).  The caller opens it around
//! any run, molecule or crystal, and closes it after; closing restores both counts.  A deck's own `threads` nests inside it.
typedef struct qcb_threads qcb_threads;
qcb_threads* qcb_threads_open(int openmp, int blas);
void         qcb_threads_close(qcb_threads* h);

//! THE SHELL WHAT-IF of a deck's basis (deck::MakeBasisShells; the GUI's conditioning table, the user 2026-11-11): the cell's
//! untrimmed basis analysed ONCE (analytic S per k, no grids), then queried per set of kept shells -- each query a few eigensolves.
typedef struct qcb_shells qcb_shells;
qcb_shells* qcb_shells_make(const char* deckJson, char* err, int errlen);
void        qcb_shells_free(qcb_shells* h);
//! {"kept": [bool per shell]} (absent: the deck's own trim) -> the library's basisTrimAnalysis RECORD (report::GlobalReport): floor,
//! healthy, lambdaMin, worstK, cond, aoDrops, nFunctions, blocks [{k, nFunctions, lambdaMin, cond, aoDrops}], shells [{Z, l,
//! exponents, nAO, kept, lambdaMinToggled, condToggled, aoDropsToggled, participation}], deckKept, autoTrimKept; or {"error"}.
//! The same two calls as qcb_materials.
int         qcb_shells_query(qcb_shells* h, const char* req, char* out, int outlen);

#ifdef __cplusplus
}
#endif
