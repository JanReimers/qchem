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

void step_trampoline(qcb_sink* s, const qcb_step* st)
{
    nb::gil_scoped_acquire gil;
    if (s->error) return;
    try { s->fn(st->iteration, st->energy, st->dE, st->commutator, st->drho); }
    catch (nb::python_error& e) { s->error.emplace(std::move(e)); }
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

} // namespace

// The Python-facing owner of one bridge handle.
class Calculator
{
public:
    Calculator(const std::vector<int>& Z, const std::vector<double>& pos, const std::string& basis,
               const std::string& model, bool symmetry, bool polarized, int multiplicity, int maxIter, const std::string& pp,
               const std::string& accelerator, int meshRadial, int meshAngular, nb::object onIteration)
    {
        if (pos.size() != 3*Z.size()) throw std::invalid_argument("positions must be flat 3N (bohr)");
        const qcb_request req{Z.data(), int(Z.size()), pos.data(), basis.c_str(), model.c_str(), pp.c_str(),
                              accelerator.c_str(), symmetry ? 1 : 0, polarized ? 1 : 0, multiplicity, maxIter,
                              meshRadial, meshAngular};
        qcb_sink sink{onIteration, std::nullopt};
        const bool stream = !onIteration.is_none();
        char err[errlen] = {0};
        {
            nb::gil_scoped_release nogil;
            h = qcb_make(&req, stream ? step_trampoline : nullptr, stream ? &sink : nullptr, err, errlen);
        }
        if (sink.error) { if (h) qcb_free(h); h = nullptr; throw std::move(*sink.error); }
        if (!h) throw std::runtime_error(err);
    }
    ~Calculator() { qcb_free(h); }
    Calculator(const Calculator&)            = delete;
    Calculator& operator=(const Calculator&) = delete;

    //! Re-run the SCF from the seed, calling on_iteration(iter, E, dE, [F,D], drho) live.  Returns converged.
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

NB_MODULE(qchem_py, mod)
{
    mod.doc() = "qchem -> Python bridge: resolve + converge a molecular SCF (qchem::Calculation), sample fields onto grids.";

    nb::class_<Calculator>(mod, "Calculator")
        .def(nb::init<const std::vector<int>&, const std::vector<double>&, const std::string&, const std::string&,
                      bool, bool, int, int, const std::string&, const std::string&, int, int, nb::object>(),
             nb::arg("numbers"), nb::arg("positions"), nb::arg("basis") = "dzvp", nb::arg("model") = "HF",
             nb::arg("symmetry") = false, nb::arg("polarized") = false, nb::arg("multiplicity") = 0, nb::arg("max_iter") = 100,
             nb::arg("pp") = "", nb::arg("accelerator") = "DIIS", nb::arg("mesh_radial") = 0, nb::arg("mesh_angular") = 0,
             nb::arg("on_iteration") = nb::none(),
             "MoleculeSetup::Resolve + converge.  Atomic numbers + flat 3N bohr positions.  on_iteration(iter, E, dE, "
             "[F,D], drho), if given, streams the SCF live (the GIL is released during the solve).")
        .def("run_scf", &Calculator::RunSCF, nb::arg("on_iteration"), nb::arg("max_iter") = 100,
             "Re-run the SCF from the seed, streaming each iteration; returns converged.")
        .def("total_energy", &Calculator::TotalEnergy)
        .def("converged",    &Calculator::Converged)
        .def("iterations",   &Calculator::Iterations)
        .def("structure",    &Calculator::Structure)
        .def("levels",       &Calculator::Levels, "Every MO, lowest first: [{'e','occ','spin','label'}].")
        .def("density",      &Calculator::Density,     nb::arg("n") = 64, nb::arg("pad") = 5.0)
        .def("spin_density", &Calculator::SpinDensity, nb::arg("n") = 64, nb::arg("pad") = 5.0,
             "m(r) = rho_up - rho_down; the zero field for an unpolarized run.")
        .def("orbital",      &Calculator::Orbital,     nb::arg("index") = 0, nb::arg("n") = 64, nb::arg("pad") = 5.0,
             "MO `index` in levels() order (0 = lowest).")
        .def("gradient", &Calculator::Gradient, nb::arg("n") = 24, nb::arg("pad") = 4.0);
}
