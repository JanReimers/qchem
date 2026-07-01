// File: pybind/qchem_py.cpp
//
// Pure nanobind layer: includes ONLY the flat C ABI (qcb_api.h) -- no qchem
// modules, so nanobind's textual std headers never clash with module GMFs.
// Allocates owned NumPy arrays and lets the C++ bridge fill them.
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>
#include "qcb_api.h"

namespace nb = nanobind;

// A NumPy array of the given shape, backed by a heap buffer it owns.
static nb::ndarray<nb::numpy, double> make_array(std::initializer_list<size_t> shape)
{
    size_t total = 1; for (size_t s : shape) total *= s;
    double* buf = new double[total];
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    return nb::ndarray<nb::numpy, double>(buf, shape, owner);
}

// Python-facing wrapper around the opaque C handle.
struct Calc { void* h = nullptr; };

NB_MODULE(qchem_py, mod)
{
    mod.doc() = "qchem -> Python bridge: run a molecular SCF, sample fields onto grids.";

    nb::class_<Calc>(mod, "Calculator")
        .def("__init__", [](Calc* self, const std::vector<int>& Z, const std::vector<double>& pos,
                            const std::string& basis, const std::string& method, int max_iter) {
                 self->h = qcb_make(Z.data(), (int)Z.size(), pos.data(),
                                    basis.c_str(), method.c_str(), max_iter);
             },
             nb::arg("numbers"), nb::arg("positions"),
             nb::arg("basis") = "dzvp", nb::arg("method") = "HF", nb::arg("max_iter") = 20,
             "Build a molecule + converge an SCF. method: HF | LDA | Xalpha.")
        .def("__del__", [](Calc& c){ if (c.h) { qcb_free(c.h); c.h = nullptr; } })
        .def("total_energy", [](Calc& c){ return qcb_energy(c.h); })
        .def("run_scf", [](Calc& c, nb::callable cb) {
                 // captureless trampoline -> C function pointer; user = &cb (alive for the call)
                 auto tramp = [](void* u, int it, double E, double dE, double comm, double drho) {
                     (*static_cast<nb::callable*>(u))(it, E, dE, comm, drho);
                 };
                 return qcb_run_scf(c.h, tramp, &cb);
             }, nb::arg("callback"),
             "Re-run the SCF from the seed, calling callback(iter, E, dE, [F,D], drho) "
             "each iteration. Returns the iteration count.")
        .def("structure", [](Calc& c) {
                 int nat = qcb_natoms(c.h);
                 std::vector<int> Z(nat); std::vector<double> xyz(3*nat);
                 qcb_atoms(c.h, Z.data(), xyz.data());
                 nb::dict d; d["numbers"] = Z; d["positions"] = xyz; return d;
             })
        .def("density", [](Calc& c, int n, double pad) {
                 auto arr = make_array({(size_t)n, (size_t)n, (size_t)n});
                 double o[3], s[3];
                 qcb_density(c.h, n, pad, arr.data(), o, s);
                 nb::dict d; d["values"] = arr;
                 d["origin"] = nb::make_tuple(o[0], o[1], o[2]);
                 d["spacing"] = nb::make_tuple(s[0], s[1], s[2]);
                 return d;
             }, nb::arg("n") = 64, nb::arg("pad") = 5.0)
        .def("orbital", [](Calc& c, int index, int n, double pad) {
                 auto arr = make_array({(size_t)n, (size_t)n, (size_t)n});
                 double o[3], s[3]; int sgn = 0;
                 qcb_orbital(c.h, index, n, pad, arr.data(), o, s, &sgn);
                 nb::dict d; d["values"] = arr;
                 d["origin"] = nb::make_tuple(o[0], o[1], o[2]);
                 d["spacing"] = nb::make_tuple(s[0], s[1], s[2]);
                 d["signed"] = (bool)sgn;
                 return d;
             }, nb::arg("index") = 0, nb::arg("n") = 64, nb::arg("pad") = 5.0)
        .def("gradient", [](Calc& c, int n, double pad) {
                 auto arr = make_array({(size_t)n, (size_t)n, (size_t)n, 3});
                 double o[3], s[3];
                 qcb_gradient(c.h, n, pad, arr.data(), o, s);
                 nb::dict d; d["vectors"] = arr;
                 d["origin"] = nb::make_tuple(o[0], o[1], o[2]);
                 d["spacing"] = nb::make_tuple(s[0], s[1], s[2]);
                 return d;
             }, nb::arg("n") = 24, nb::arg("pad") = 4.0);
}
