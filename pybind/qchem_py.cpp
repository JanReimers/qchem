// File: pybind/qchem_py.cpp
//
// Pure nanobind layer: includes ONLY the flat C ABI (qcb_api.h) -- no qchem modules, so nanobind's textual std
// headers never clash with module GMFs.  Allocates owned NumPy arrays and lets the C++ bridge fill them.
//
// THREADING: the SCF runs with the GIL RELEASED (so a GUI worker thread does not freeze the window); the step
// trampoline re-acquires it only to call the Python callable.  A Python exception raised by that callable is
// held in the sink, later steps are skipped, and it is re-raised once the SCF returns.
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>
#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#include "qcb_api.h"

namespace nb = nanobind;

// The sink the bridge streams to: the Python callable, plus the first error it raised.
struct qcb_sink
{
    nb::object                      fn;
    std::optional<nb::python_error> error;
};

namespace {

//! on_iteration(iteration, energy, cells): cells is a list of (title, value, text, tolerance | None, "log" | "linear"),
//! one per console column in display order -- generic, so a column the lib adds needs no change here.  An exception the
//! callable raises STOPS the SCF (the GUI's Stop button raises one) and is re-raised to the caller once the solve unwinds.
int step_trampoline(qcb_sink* s, const qcb_row* r)
{
    nb::gil_scoped_acquire gil;
    if (s->error) return 1;
    try
    {
        nb::list cells;
        for (int k = 0; k < r->ncells; ++k)
        {
            const qcb_cell& c = r->cells[k];
            cells.append(nb::make_tuple(nb::str(c.title), c.value, nb::str(c.text),
                                        c.tolerance == c.tolerance ? nb::object(nb::float_(c.tolerance)) : nb::none(),
                                        c.logScale ? "log" : "linear"));
        }
        s->fn(r->iteration, r->energy, cells);
        return 0;
    }
    catch (nb::python_error& e) { s->error.emplace(std::move(e)); }
    return 1;
}

// A NumPy array of the given shape, backed by a heap buffer it owns.
nb::ndarray<nb::numpy, double> make_array(std::initializer_list<size_t> shape)
{
    size_t total = 1; for (size_t s : shape) total *= s;
    double* buf = new double[total];
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    return nb::ndarray<nb::numpy, double>(buf, shape, owner);
}

nb::dict grid_dict(const char* key, nb::ndarray<nb::numpy, double> arr, const double* o, const double* s)
{
    nb::dict d; d[key] = arr;
    d["origin"]  = nb::make_tuple(o[0], o[1], o[2]);
    d["spacing"] = nb::make_tuple(s[0], s[1], s[2]);
    return d;
}

constexpr int errlen = 1024;

//! A saved state the facade REFUSED to restart from (kind, basis, electrons, blocks, or an unreadable file): a value the caller acts
//! on (a fresh run, or another state), raised as qchem_py.RestartRefused (a RuntimeError).
struct RestartRefused : std::runtime_error { using std::runtime_error::runtime_error; };

//! A bridge call that writes JSON text (length-returning, sized by a first call with outlen 0) -> the Python object.
template <class F> nb::object json_call(F&& f)
{
    std::string out(size_t(f(nullptr, 0)) + 1, '\0');
    out.resize(size_t(f(out.data(), int(out.size()))));
    return nb::module_::import_("json").attr("loads")(nb::str(out.c_str(), out.size()));
}

} // namespace

// The Python-facing owner of one bridge handle.
class Calculator
{
public:
    Calculator(const std::vector<int>& Z, const std::vector<double>& pos, const std::string& basis,
               const std::string& model, bool symmetry, bool polarized, int multiplicity, int maxIter, const std::string& pp,
               const std::string& accelerator, int meshRadial, int meshAngular, double relaxStart, double relaxMax,
               int diisDepth, double diisSVTol, double diisEngageFD, bool mom, int momStart, double relaxMin,
               double ladderSwitchAt, double gateFD, double gateDrho, double gateDE, double gateDFD,
               const std::string& mixer, double smearingkT, const std::string& restartFrom, bool buildOnly, nb::object onIteration)
    {
        if (pos.size() != 3*Z.size()) throw std::invalid_argument("positions must be flat 3N (bohr)");
        const qcb_request req{Z.data(), int(Z.size()), pos.data(), basis.c_str(), model.c_str(), pp.c_str(),
                              accelerator.c_str(), symmetry ? 1 : 0, polarized ? 1 : 0, multiplicity, maxIter,
                              meshRadial, meshAngular, relaxStart, relaxMax, diisDepth, diisSVTol, diisEngageFD,
                              mom ? 1 : 0, momStart, relaxMin, ladderSwitchAt,
                              gateFD, gateDrho, gateDE, gateDFD, mixer.c_str(), smearingkT};
        qcb_sink sink{onIteration, std::nullopt};
        const bool stream = !onIteration.is_none();
        char err[errlen] = {0};
        int refused = 0;
        {
            nb::gil_scoped_release nogil;
            h = qcb_start(&req, restartFrom.c_str(), buildOnly ? 1 : 0, stream ? step_trampoline : nullptr, stream ? &sink : nullptr, &refused, err, errlen);
        }
        if (sink.error) { if (h) qcb_free(h); h = nullptr; throw std::move(*sink.error); }
        if (!h && refused) throw RestartRefused(err);
        if (!h) throw std::runtime_error(err);
    }
    ~Calculator() { qcb_free(h); }
    Calculator(const Calculator&)            = delete;
    Calculator& operator=(const Calculator&) = delete;

    //! Re-run the SCF from the seed, calling on_iteration(iter, E, cells) live.  Returns converged (failure() says why not).
    bool RunSCF(nb::object onIteration, int maxIter)
    {
        qcb_sink sink{onIteration, std::nullopt};
        char err[errlen] = {0};
        int rc;
        {
            nb::gil_scoped_release nogil;
            rc = qcb_converge(h, maxIter, step_trampoline, &sink, err, errlen);
        }
        if (sink.error) throw std::move(*sink.error);
        if (rc < 0) throw std::runtime_error(err);
        return rc == 1;
    }

    //! Why the last SCF did not converge (the facade's SCFFailure as one line); None when it converged.
    nb::object Failure() const
    {
        char out[errlen] = {0};
        return qcb_failure(h, out, errlen) ? nb::object(nb::str(out)) : nb::none();
    }

    //! SAVE the run's state (HDF5 "qchem-state 2.0"), converged or not.
    void SaveState(const std::string& path) const
    {
        char err[errlen] = {0};
        int ok;
        { nb::gil_scoped_release nogil; ok = qcb_save_state(h, path.c_str(), err, errlen); }
        if (!ok) throw std::runtime_error(err);
    }
    nb::object JudgeState(const std::string& path) const
    {
        return json_call([&](char* out, int n){ nb::gil_scoped_release nogil; return qcb_judge_state(h, path.c_str(), out, n); });
    }
    nb::object RestartedFrom() const { return json_call([&](char* out, int n){ return qcb_restarted_from(h, out, n); }); }

    double TotalEnergy() const { return qcb_energy(h); }
    bool   Converged()   const { return qcb_converged(h) == 1; }
    int    Iterations()  const { return qcb_iterations(h); }

    nb::dict Structure() const
    {
        const int nat = qcb_natoms(h);
        std::vector<int> Z(nat); std::vector<double> xyz(3*nat);
        qcb_atoms(h, Z.data(), xyz.data());
        nb::dict d; d["numbers"] = Z; d["positions"] = xyz; return d;
    }

    //! Every MO, lowest first: a list of {"e", "occ", "spin" (-1/0/+1), "label"}.
    nb::list Levels() const
    {
        constexpr int nlabel = 64;
        const int n = qcb_num_levels(h);
        std::vector<double> e(n), occ(n); std::vector<int> spin(n); std::vector<char> label(size_t(n)*nlabel, '\0');
        qcb_levels(h, e.data(), occ.data(), spin.data(), label.data(), nlabel);
        nb::list out;
        for (int i = 0; i < n; ++i)
        {
            nb::dict d; d["e"] = e[i]; d["occ"] = occ[i]; d["spin"] = spin[i];
            d["label"] = std::string(label.data() + size_t(i)*nlabel);
            out.append(d);
        }
        return out;
    }

    nb::dict Density(int n, double pad) const
    {
        auto arr = make_array({size_t(n), size_t(n), size_t(n)});
        double o[3], s[3];
        qcb_density(h, n, pad, arr.data(), o, s);
        return grid_dict("values", arr, o, s);
    }

    //! The current iterate's density (spin: its m(r)); see qcb_iterate_density for WHEN it may be called.
    nb::dict IterateDensity(int n, double pad, bool spin) const
    {
        auto arr = make_array({size_t(n), size_t(n), size_t(n)});
        double o[3], s[3];
        char err[errlen] = {0};
        int ok;
        { nb::gil_scoped_release nogil; ok = qcb_iterate_density(h, spin ? 1 : 0, n, pad, arr.data(), o, s, err, errlen); }
        if (!ok) throw std::runtime_error(err);
        return grid_dict("values", arr, o, s);
    }

    nb::dict SpinDensity(int n, double pad) const
    {
        auto arr = make_array({size_t(n), size_t(n), size_t(n)});
        double o[3], s[3];
        qcb_spin_density(h, n, pad, arr.data(), o, s);
        return grid_dict("values", arr, o, s);
    }

    nb::dict Orbital(int index, int n, double pad) const
    {
        auto arr = make_array({size_t(n), size_t(n), size_t(n)});
        double o[3], s[3];
        if (!qcb_orbital(h, index, n, pad, arr.data(), o, s))
            throw std::out_of_range("orbital index " + std::to_string(index) + " is not a level");
        nb::dict d = grid_dict("values", arr, o, s);
        d["signed"] = true;
        return d;
    }

    nb::dict Gradient(int n, double pad) const
    {
        auto arr = make_array({size_t(n), size_t(n), size_t(n), 3});
        double o[3], s[3];
        qcb_gradient(h, n, pad, arr.data(), o, s);
        return grid_dict("vectors", arr, o, s);
    }

private:
    qcb_calc* h = nullptr;
};

// The Python-facing owner of one periodic calculation, built from a deck (qchem::deck::Build).
class SolidRun
{
public:
    explicit SolidRun(const std::string& deckJson)
    {
        char err[errlen] = {0};
        {
            nb::gil_scoped_release nogil;
            h = qcb_solid_make(deckJson.c_str(), err, errlen);
        }
        if (!h) throw std::runtime_error(err);
    }
    ~SolidRun() { qcb_solid_free(h); }
    SolidRun(const SolidRun&)            = delete;
    SolidRun& operator=(const SolidRun&) = delete;

    //! Converge the deck's schedule, calling on_iteration(iter, E, cells) live.  Returns converged (failure() says why not).
    bool Converge(nb::object onIteration)
    {
        qcb_sink sink{onIteration, std::nullopt};
        const bool stream = !onIteration.is_none();
        char err[errlen] = {0};
        int rc;
        {
            nb::gil_scoped_release nogil;
            rc = qcb_solid_converge(h, stream ? step_trampoline : nullptr, stream ? &sink : nullptr, err, errlen);
        }
        if (sink.error) throw std::move(*sink.error);
        if (rc < 0) throw std::runtime_error(err);
        return rc == 1;
    }

    nb::object Failure() const
    {
        char out[errlen] = {0};
        return qcb_solid_failure(h, out, errlen) ? nb::object(nb::str(out)) : nb::none();
    }
    //! None unless the last SCF converged.
    nb::object TotalEnergy() const
    {
        const double e = qcb_solid_energy(h);
        return e == e ? nb::object(nb::float_(e)) : nb::none();
    }
    int Iterations() const { return qcb_solid_iterations(h); }
    nb::object JudgeState(const std::string& path) const
    {
        return json_call([&](char* out, int n){ nb::gil_scoped_release nogil; return qcb_solid_judge_state(h, path.c_str(), out, n); });
    }
    nb::object RestartedFrom() const { return json_call([&](char* out, int n){ return qcb_solid_restarted_from(h, out, n); }); }
    std::string Resolved() const
    {
        std::string out(size_t(qcb_solid_resolved(h, nullptr, 0)) + 1, '\0');
        out.resize(size_t(qcb_solid_resolved(h, out.data(), int(out.size()))));
        return out;
    }

private:
    qcb_solid* h = nullptr;
};

NB_MODULE(qchem_py, mod)
{
    mod.doc() = "qchem -> Python bridge: a molecular SCF (qchem::MoleculeCalculation, fields sampled onto grids) and a periodic one "
                "built from an input deck (qchem::deck::Build).";

    nb::exception<RestartRefused>(mod, "RestartRefused", PyExc_RuntimeError);

    nb::class_<Calculator>(mod, "Calculator")
        .def(nb::init<const std::vector<int>&, const std::vector<double>&, const std::string&, const std::string&,
                      bool, bool, int, int, const std::string&, const std::string&, int, int, double, double,
                      int, double, double, bool, int, double, double, double, double, double, double,
                      const std::string&, double, const std::string&, bool, nb::object>(),
             nb::arg("numbers"), nb::arg("positions"), nb::arg("basis") = "dzvp", nb::arg("model") = "HF",
             nb::arg("symmetry") = false, nb::arg("polarized") = false, nb::arg("multiplicity") = 0, nb::arg("max_iter") = 100,
             nb::arg("pp") = "", nb::arg("accelerator") = "DIIS", nb::arg("mesh_radial") = 0, nb::arg("mesh_angular") = 0,
             nb::arg("relax_start") = 0.0, nb::arg("relax_max") = 0.0, nb::arg("diis_depth") = 0, nb::arg("diis_sv_tol") = 0.0,
             nb::arg("diis_engage_fd") = 0.0, nb::arg("mom") = false, nb::arg("mom_start") = 0,
             nb::arg("relax_min") = 0.0, nb::arg("ladder_switch_at") = -1.0,
             nb::arg("gate_fd") = -1.0, nb::arg("gate_drho") = -1.0, nb::arg("gate_de") = -1.0, nb::arg("gate_dfd") = -1.0,
             nb::arg("mixer") = "Auto", nb::arg("smearing_kT") = 0.0, nb::arg("restart_from") = "", nb::arg("build_only") = false,
             nb::arg("on_iteration") = nb::none(),
             "MoleculeSetup::Resolve + converge.  Atomic numbers + flat 3N bohr positions.  on_iteration(iter, E, cells), "
             "if given, streams the SCF live (the GIL is released during the solve); cells = [(title, value, text, "
             "tolerance|None, 'log'|'linear')], the console's columns in display order; an exception it raises stops the SCF and "
             "is re-raised.  relax_start / relax_max: linear D-mixing's starting step and ceiling; diis_depth / diis_sv_tol / "
             "diis_engage_fd: DIIS's history, singular-value floor and engage [F,D] (also the Ladder's DIIS rung); mom / "
             "mom_start: the Maximum Overlap Method and its reference-capture iteration (0 = the lib's default throughout).  relax_min: the "
             "adaptive linear step's floor.  ladder_switch_at: the Ladder's DIIS -> GDM hand-off at |dE/E| below it (0 = stall only, "
             "< 0 = the lib's 1e-6).  gate_fd / gate_drho / gate_de / gate_dfd: the convergence gates [F,D], Δρ per electron, "
             "|ΔE/E|, |Δ[F,D]| (> 0 the gate, 0 = off, < 0 = the lib's default).  mixer: SCFParams::Mixing by name (Auto | Linear | "
             "Kerker | Pulay; the facade refuses what a molecule cannot mix).  smearing_kT: Fermi-Dirac smearing, Ha (0 = off).  "
             "restart_from: a saved state (save_state's file) to start from (MoleculeCalculation::Restart): an exact resume or a "
             "warm start (restarted_from() says which); a state the facade refuses raises RestartRefused.  build_only: BUILD the whole graph (basis, every mesh, the SCF "
             "iterator) and stop before iteration 1, filed as the report run '<ID> build'; run_scf then converges on that graph.")
        .def("run_scf", &Calculator::RunSCF, nb::arg("on_iteration"), nb::arg("max_iter") = 100,
             "Run the SCF (after build_only: on the built graph; else again from the seed), streaming each iteration; returns converged.")
        .def("save_state", &Calculator::SaveState, nb::arg("path"),
             "Save the run's state (HDF5 'qchem-state 2.0'), converged or not; a later Calculator(..., restart_from=path) starts from it.")
        .def("judge_state", &Calculator::JudgeState, nb::arg("path"),
             "What a restart of this run from `path` would do, without doing it: {'exact', 'warm': [line, ...]} or "
             "{'refused': {'why', 'details'}}.")
        .def("restarted_from", &Calculator::RestartedFrom,
             "What this run's restart did: {'path', 'exact', 'warm', 'maxDisplacement' (bohr), 'electrons'}; None for a fresh run.")
        .def("total_energy", &Calculator::TotalEnergy)
        .def("converged",    &Calculator::Converged)
        .def("failure",      &Calculator::Failure, "Why the last SCF did not converge (SCFFailure: why, iterations, "
             "details); None when it converged.")
        .def("iterations",   &Calculator::Iterations)
        .def("structure",    &Calculator::Structure)
        .def("levels",       &Calculator::Levels, "Every MO, lowest first: [{'e','occ','spin','label'}].")
        .def("density",      &Calculator::Density,     nb::arg("n") = 64, nb::arg("pad") = 5.0)
        .def("spin_density", &Calculator::SpinDensity, nb::arg("n") = 64, nb::arg("pad") = 5.0,
             "m(r) = rho_up - rho_down; the zero field for an unpolarized run.")
        .def("orbital",      &Calculator::Orbital,     nb::arg("index") = 0, nb::arg("n") = 64, nb::arg("pad") = 5.0,
             "MO `index` in levels() order (0 = lowest).")
        .def("iterate_density", &Calculator::IterateDensity, nb::arg("n") = 64, nb::arg("pad") = 5.0, nb::arg("spin") = false,
             "The CURRENT SCF iterate's density (spin=True: its spin density), sampled like density(): the latest diagonalisation's, "
             "built fresh per call.  Call it ONLY from run_scf's on_iteration callback (the compute thread, between iterations) or "
             "between runs; raises before the first iteration.")
        .def("gradient", &Calculator::Gradient, nb::arg("n") = 24, nb::arg("pad") = 4.0);

    mod.def("report", []()
    {
        std::string out(size_t(qcb_report(nullptr, 0)) + 1, '\0');
        out.resize(size_t(qcb_report(out.data(), int(out.size()))));
        return out;
    }, "The library's run record (JSON text): every run since the last clear_report, keyed \"<name>@<start>\".");
    mod.def("current_report", []()
    {
        std::string out(size_t(qcb_current_report(nullptr, 0)) + 1, '\0');
        out.resize(size_t(qcb_current_report(out.data(), int(out.size()))));
        return out;
    }, "The open run's record so far (JSON text, \"\" outside a run): call it from on_iteration, on the run's own thread.");
    mod.def("clear_report", &qcb_clear_report, "Drop the accumulated run record (a no-op while a run is open).");
    mod.def("point_group", [](const std::vector<int>& species, const std::vector<double>& pos, double tol)
    {
        if (pos.size() != 3*species.size()) throw std::invalid_argument("positions must be flat 3N (bohr)");
        char sym[64] = {0}, ab[64] = {0};
        const int order = qcb_point_group(species.data(), int(species.size()), pos.data(), tol, sym, 64, ab, 64);
        if (order < 0) throw std::runtime_error(sym);
        nb::dict d; d["symbol"] = std::string(sym); d["abelian"] = std::string(ab); d["order"] = order; return d;
    }, nb::arg("species"), nb::arg("positions"), nb::arg("tol") = 1e-3,
       "The DETECTED molecular point group (Symmetry::Molecule::DetectPointGroup): {'symbol' (Schoenflies), 'abelian' "
       "(the SCF's blocking subgroup), 'order' (0 for the linear groups)}.  species: one integer label per atom (Z, or Z "
       "combined with a moment sign to see a magnetic decoration).");
    mod.def("space_group", [](const std::vector<double>& lattice, const std::vector<int>& species,
                              const std::vector<double>& frac, double tol)
    {
        if (lattice.size() != 9 || frac.size() != 3*species.size())
            throw std::invalid_argument("lattice must be 9 (rows a, b, c) and frac flat 3N");
        constexpr int maxops = 192;
        std::vector<double> W(9*maxops);
        int symmorphic = 0;
        char err[errlen] = {0};
        const int n = qcb_space_group(lattice.data(), species.data(), int(species.size()), frac.data(), tol,
                                      W.data(), maxops, &symmorphic, err, errlen);
        if (n < 0) throw std::runtime_error(err);
        W.resize(size_t(9*std::min(n, maxops)));
        nb::dict d; d["order"] = n; d["symmorphic"] = symmorphic == 1; d["W"] = W; return d;
    }, nb::arg("lattice"), nb::arg("species"), nb::arg("frac"), nb::arg("tol") = 1e-4,
       "The DETECTED space group (Symmetry::Lattice_3D::SpaceGroup::Detect) of a primitive cell: {'order', 'symmorphic', "
       "'W' (each op's linear part in lattice coordinates, flat row-major 9 per op)}.");

    mod.def("materials", [](const std::string& req)
    {
        int n;
        { nb::gil_scoped_release nogil; n = qcb_materials(req.c_str(), nullptr, 0); }
        std::string out(size_t(n) + 1, '\0');
        out.resize(size_t(qcb_materials(nullptr, out.data(), int(out.size()))));
        return out;
    }, nb::arg("request"),
       "Crystal files and the materials layers (JSON text in, JSON text out; {\"error\": ...} on a refusal): list | import | "
       "export | save | cell -- see qcb_api.h qcb_materials (the caller passes the user / workspace library paths).");

    nb::class_<SolidRun>(mod, "SolidRun")
        .def(nb::init<const std::string&>(), nb::arg("deck_json"),
             "deck::Build: parse the input deck (JSON text, schema 1 -- the file rundeck runs) and build its periodic "
             "calculation, stopping before the SCF.  Raises RuntimeError with the library's text on a refusal (an unknown "
             "key, an unknown structure, a setup the library refuses).")
        .def("converge", &SolidRun::Converge, nb::arg("on_iteration") = nb::none(),
             "Converge the deck's schedule (deck::ScheduleToRun), on_iteration(iter, E, cells) live as for Calculator; "
             "returns converged.")
        .def("failure",      &SolidRun::Failure, "Why the last SCF did not converge; None when it converged.")
        .def("total_energy", &SolidRun::TotalEnergy, "The converged total energy (hartree); None unless converged.")
        .def("iterations",   &SolidRun::Iterations)
        .def("judge_state",  &SolidRun::JudgeState, nb::arg("path"), "As Calculator.judge_state, for this crystal.")
        .def("restarted_from", &SolidRun::RestartedFrom,
             "What the deck's state.restartFrom did (as Calculator.restarted_from); None for a fresh run.")
        .def("resolved",     &SolidRun::Resolved, "The resolved deck (JSON text): every default filled in, what rundeck records.");
}
