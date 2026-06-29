"""
qviz.backend_qchem -- the REAL ComputeBackend, backed by the nanobind module.

This replaces backend_analytic.AnalyticBackend with genuine qchem output: it
runs a molecular HF SCF in C++ and samples the converged density / HOMO /
gradient onto grids. The frontends (PyVista, the desktop app) are unchanged --
they still only see qviz.data records.

Requires the compiled extension `qchem_py` (built from pybind/ with
-DQCHEM_PYBIND=ON). We add the PIC build dir to sys.path so no install step is
needed during development.
"""
from __future__ import annotations
import sys, pathlib
import numpy as np

from .data import Structure, ScalarField, VectorField, SCFStep, ComputeBackend

# repo root = .../qchem6 ; the module lands in build/PIC/pybind/
_REPO = pathlib.Path(__file__).resolve().parents[2]
for _cand in (_REPO / "build/PIC/pybind", _REPO / "build/Release/pybind"):
    if _cand.is_dir():
        sys.path.insert(0, str(_cand))
        break

import qchem_py   # the nanobind extension

# minimal Z -> symbol (extend as needed; the C bridge returns Z, symbols map here)
_SYMBOL = {1: "H", 3: "Li", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na",
           14: "Si", 15: "P", 16: "S", 17: "Cl", 25: "Mn", 27: "Co", 28: "Ni"}


# water, experimental geometry in BOHR (matches UnitTests/M_HF_U.C MakeWater)
WATER_NUMBERS = [8, 1, 1]
WATER_POSITIONS = [0.0, 0.0, 0.0,
                   0.0, 1.431, 1.107,
                   0.0, -1.431, 1.107]


class QChemBackend(ComputeBackend):
    def __init__(self, numbers=None, positions=None, basis="dzvp", max_iter=20):
        numbers = list(numbers if numbers is not None else WATER_NUMBERS)
        positions = list(np.asarray(positions if positions is not None
                                    else WATER_POSITIONS, float).ravel())
        self._calc = qchem_py.Calculator(numbers, positions, basis, max_iter)
        self._struct = self._read_structure()

    # -- geometry ----------------------------------------------------------
    def _read_structure(self) -> Structure:
        d = self._calc.structure()
        numbers = np.asarray(d["numbers"], int)
        pos = np.asarray(d["positions"], float).reshape(-1, 3)
        return Structure(symbols=[_SYMBOL.get(int(z), "X") for z in numbers],
                         positions=pos, numbers=numbers)

    def structure(self) -> Structure:
        return self._struct

    def total_energy(self) -> float:
        return self._calc.total_energy()

    # -- fields ------------------------------------------------------------
    @staticmethod
    def _scalar(d, name, signed) -> ScalarField:
        return ScalarField(name, np.asarray(d["values"]),
                           tuple(d["origin"]), tuple(d["spacing"]), signed=signed)

    def density(self, n: int = 64, pad: float = 5.0) -> ScalarField:
        return self._scalar(self._calc.density(n, pad), "Electron density", False)

    def orbital(self, index: int = 0, n: int = 64, pad: float = 5.0) -> ScalarField:
        d = self._calc.orbital(index, n, pad)
        return self._scalar(d, f"MO #{index} (HOMO-{index})" if index else "HOMO", bool(d["signed"]))

    def density_gradient(self, n: int = 22, pad: float = 4.0) -> VectorField:
        d = self._calc.gradient(n, pad)
        return VectorField("grad(rho)", np.asarray(d["vectors"]),
                           tuple(d["origin"]), tuple(d["spacing"]))

    # -- SCF: static-first stub --------------------------------------------
    def run_scf(self):
        """Live streaming is the next pass (needs an observer hook in
        SCFIterator). For now yield a single converged point so the desktop
        app's convergence panel shows the final energy without erroring."""
        yield SCFStep(iteration=0, energy=self.total_energy(),
                      dE=0.0, commutator=0.0, drho=0.0)
