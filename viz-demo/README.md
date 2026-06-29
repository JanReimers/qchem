# qchem viz — feasibility demo

A spike for the qchem6 GUI layer. It proves out the recommended stack and, more
importantly, the **seam** that keeps the C++ science decoupled from the plots:

```
   qchem (C++)  ──nanobind──▶  qviz.ComputeBackend  ──▶  PyVista (3D) / pyqtgraph (2D)
        │                            │                          │
   already exists          the only contract             flashy frontends
                           the GUI knows about
```

Today the backend is `qviz.backend_analytic.AnalyticBackend` (closed-form water:
density, a HOMO, gradient, a simulated SCF). Swap in a `backend_qchem.QChemBackend`
nanobind module that returns the *same records* from `src/`, and every frontend
keeps working untouched.

## The seam — `qviz/data.py`

Frontend-agnostic records + a `ComputeBackend` Protocol:

| record          | feeds                          | nanobind source in qchem            |
|-----------------|--------------------------------|-------------------------------------|
| `Structure`     | ball-and-stick                 | Molecule / Structure geometry       |
| `ScalarField`   | isosurfaces, slices, volume    | ChargeDensity / orbital on a grid   |
| `VectorField`   | glyphs (∇ρ, forces, current)   | density gradient / forces           |
| `SCFStep`       | live convergence plot          | yielded from the SCF iterate loop   |

Arrays are contiguous NumPy → zero-copy views of Blaze storage across nanobind.
That zero-copy is what makes per-iteration live plotting cheap.

## Run

```bash
uv venv -p 3.12 .venv && . .venv/bin/activate   # Python 3.12: VTK/PySide wheels exist
uv pip install -r requirements.txt

python scripts/make_artifacts.py    # headless: renders everything to out/
python app_desktop.py               # the integrated desktop shell (needs a display)
```

**Wayland note:** VTK's interactor needs X11/GLX, which a native Wayland surface
can't host (`BadWindow` on `X_ConfigureWindow`). `app_desktop.py` auto-forces
`QT_QPA_PLATFORM=xcb` (XWayland) when it detects a Wayland session — no action
needed. To override: `QT_QPA_PLATFORM=wayland python app_desktop.py`.

## What's here

- `qviz/data.py` — the C++/Python boundary (records + `ComputeBackend`).
- `qviz/backend_analytic.py` — **REPLACE-ME** stand-in for `nanobind(qchem)`.
- `qviz/scene.py` — the only PyVista-aware module (records → actors).
- `qviz/project.py` — workspace save/restore on HDF5 (read in C++ too, via HighFive).
- `app_desktop.py` — PySide6 + pyvistaqt + pyqtgraph shell: 3D viewport, live SCF
  plot, iso-level + slice sliders, project open/save.
- `scripts/make_artifacts.py` — off-screen renders → `out/` (works without a GPU/display\*).

\* SSAA anti-aliasing is skipped automatically on GPU-less boxes.

## out/ artifacts

`density_iso.png`, `orbital_homo.png`, `gradient_field.png`,
`slice_sweep.gif` (slider sweeping a 2D slice through the 3D density),
`scf_convergence.gif` (streamed ΔE / ‖[F,D]‖ / ‖Δρ‖), `contact_sheet.png`,
and `water.qproj.h5` (a saved project, round-tripped on load).

## Next step toward the real thing

Write `backend_qchem.QChemBackend(ComputeBackend)` with nanobind:
1. expose `Molecule` geometry → `Structure`;
2. sample an existing `ChargeDensity` on a grid → `ScalarField` (zero-copy);
3. add a per-iteration callback to the SCF driver that `yield`s `SCFStep`.
Point `app_desktop.py` at it instead of `AnalyticBackend`. Done.
