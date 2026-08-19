# The qchem ↔ CP2K head-to-head — runtime, RAM, energy

**Step 1 of `doc/OpenWork.md`.  This is an INSTRUMENT, not a report.**  Steps 2–4 (folds, Φ-screening, RAM)
have no finish line without it: everything optimised in the 2026-08 runtime campaign was measured against a
single hand-run disabled test at `GPW_MNO_NMAX=3`, which is how two charter premises survived for months
while being wrong.  A row here is only meaningful with its PROVENANCE columns — same span held full-rank by
both codes, and the fold/threading/BLAS state that produced the number.

## How to produce a qchem row

```bash
GPW_REPORT=1 build/Release/IntegrationTests/ITMain --gtest_filter=<test> --gtest_also_run_disabled_tests
```

Every GPW run now prints, without extra flags:

| what | where it comes from |
|---|---|
| `Etot`, iterations, convergence verdict | the run fingerprint line |
| wall time, and the per-bucket ledger | the `timing` section |
| **PEAK RSS (MB, process high-water)** | `timing`, from Linux `VmHWM` (added 2026-08-19) |
| `[fold] <site>: … = F×` for every fold site | `EmitFold`, Step 0b |
| `[stream cache] … pts64/pts32/runs/meanRun/dropped` | the stream-cache readout |
| `[site moments] … [e]` (polarized) | `QCHEM_SITE_MOMENTS=1`, Step 0a |

**Read the provenance, not just the number.**  `PEAK RSS` is the process high-water mark, so it is only a
clean per-config figure when the process runs ONE config — a `--gtest_filter` naming several tests reports
the watermark of the whole process.  Thread count is `GPW_OMP_THREADS` (serial when unset); BLAS routing is
`QCHEM_BLAZE_BLAS` (default ON).

## The table

Energies in Ha.  qchem measured on this box (14 GB, 16 cores); CP2K numbers are the banked oracles from
`doc/CP2Kresults.md` (CP2K 2026.1 serial `ssmp`, built at `~/Code/cp2k` — **not installed on the dev box, so
the CP2K column is filled from banked runs, never measured here**).

| system | k-mesh | span | qchem Etot | CP2K Etot | Δ | qchem wall | CP2K wall | qchem peak RSS | CP2K peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| Si (FCC) | Γ | SIPP_SR | −7.11507 | −7.11506 | **1e-5** | 6.1 s | 3.5 s | 267 MB | ❓ |
| Si (FCC) | 2×2×2 Γ-centred | SIPP_SR | −7.8 ¹ | −7.86744 | ¹ | 8.5 s | 4.8 s | 269 MB | ❓ |
| NaF (rocksalt) | Γ | VALENCE-LOWQ-SR | ❓ | −27.93128 ² | | ❓ | 9m32s | ❓ | ❓ |
| MnO AFM-II | Γ | VA (N=118) | −61.40298 ³ | −61.30333 | **−99.7 mHa** | ❓ | minutes ³ | ❓ | ❓ |
| MnO FM | Γ | VA (N=118) | −61.44158 ³ | −61.30478 | **−136.8 mHa** | ❓ | minutes ³ | ❓ | ❓ |
| MnO AFM-II | 2×2×2 (`MNO_KMESH=2`) | VA | ❓ | ❓ | | ❓ | ❓ | ❓ | ❓ |

¹ the fingerprint prints 2 s.f. at default detail; re-run for the full digits before quoting a Δ.
² the −27.93128 "oracle" was RETRACTED as a 3.5 Ha `EPS_PGF_ORB` screening artifact — the tight-eps rerun
gives −24.432.  See `doc/GPWPlan1.md`; the NaF row needs re-basing on the tight-eps number.
³ from the VA verdict table in `doc/SphericalLatticePlan.md` (runs 61/62), which recorded ENERGIES only.

### What is missing, and who has to produce it

- **Every CP2K RAM figure.**  CP2K does not print peak RSS by default; the runs must be wrapped
  (`/usr/bin/time -v`, or `systemd-run --scope` + `MemoryPeak`).  Without this the RAM half of Steps 2/4 has
  no reference at all.
- **MnO wall/RAM on both sides** — the VA runs banked energies only.
- **NaF** on both sides, re-based on the tight-eps CP2K number.
- **The multi-k MnO row** (`MNO_KMESH=2`, ALL-TRIM so the whole k-sum runs real) — cost is unmeasured.

### Like-for-like is available today

CP2K can be forced down to the spherical spans qchem holds at FULL RANK — the decks and basis entries
already exist: `IntegrationTests/CP2K/mno_{afm2,fm}_gpw_v{a,b}.inp` with the `VALENCE-LOWQ-V{A,B}` entries in
`IntegrationTests/CP2K/VALENCE-LOWQ-BASIS`; **VA = N 118**, **VB = N 128**, and qchem has run VA at full rank
("kept 118 of 118").  So the MnO comparison does NOT wait on the 136-span question (`OpenWork` Step 6).
Si and NaF need the same treatment: **one shared span per material, held full-rank by both codes, or the row
is not a comparison.**

## Standing observations

- **Si Γ agrees to 1e-5 Ha** and costs ~1.7× CP2K's wall time at 267 MB.  On the small cells the runtime gap
  is modest; it is the diffuse/magnetic cells where it opens up.
- **MnO carries a ~100 mHa configuration-BLIND offset and a ~37 mHa configuration-SELECTIVE one** on an
  identical span (`OpenWork` Step 5) — and qchem sits BELOW a variational reference there, which convicts an
  operator or convention rather than a basis.
- **The production MnO run folds NOTHING** — `[fold]` prints `NONE` at all three sites on a cell whose
  magnetic group has 12–24 ops (`OpenWork` Step 2).  Any MnO runtime row taken before that changes is
  measuring an unfolded run, and must say so.
