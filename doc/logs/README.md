# doc/logs — full run logs, kept for reading

**Save the WHOLE log here, never a tail.**  The user reads these and has caught real findings in them
(user, 2026-08-25).  A truncated capture also silently discards ctest's `100% tests passed …` summary and,
if you pipe through `tail`, reports the PIPELINE's exit code instead of the tool's — both bit this session.

`scripts/bench <label> -- <cmd>` writes `bench_<slug>.log` plus `bench_<slug>.log.time`
(`/usr/bin/time -v`: peak RSS, CPU%, wall).  Move both here when the run is worth keeping.

`old/` — everything from before 2026-08-25 (the MnO campaign runs 1–65 and the 2026-08-19/20 benchmark
cuts).  Kept, not curated.

## 2026-08-25 — the `GPW_XC_DM_SOURCE` promotion test (doc/OpenWork.md item 1)

All rows: MnO AFM-II, VA (N=118), `MNO_IMPOSE=1`, **single thread**, one binary (the banked-row build).
Full commands are in each log's header.

| log | arm | Etot | Eee | m_stag |
|---|---|---|---|---|
| `st_baseline` | default path — **the control** | **−61.40297621** | 13.480 | survived |
| `st_rankcensus` | default + `GPW_DM_RANK=1` (pivot spectrum vs kT) | −61.40297621 | — | survived |
| `st_siteblockfix` | default, AFTER the site-block fix + `QCHEM_SITE_MOMENTS=1` | −61.40297621 | 13.480 | survived; **Mn ±4.78 e** |
| `st_pulay8_control` | default + `MNO_PULAY=8` | −61.40297621 | 13.480 | survived |
| `st_dmsource` | `GPW_XC_DM_SOURCE=1` | −45.52875429 | **28.995** | died it 12 |
| `st_dmsource_boost2` | + `GPW_XC_DM_BOOST=2` (killed in stage 2) | −45.36 (st1) | — | died it 10 |
| `st_dmsource_boost05` | + `GPW_XC_DM_BOOST=0.5` | −45.52847324 | — | died it 19 |
| `st_dmsource_fullrank` | + `QCHEM_DM_LOWRANK=0` | −45.52875429 | — | died it 12 |
| `st_dmsource_pulay8` | + `MNO_PULAY=8` | −46.32003443 | 30.060 | st1 survived, st2 collapsed |
| `st_flatkerker` | **no flag**, `MNO_KERKER_G0=0.01` | **−38.45368688** | **35.098** | **died it 7** |
| `st_nokerker` | **no flag**, `MNO_KERKER_G0=0` (→ linear D-mix) | −56.38756586 | 13.609 | died it 13 |
| `st_held_armed` / `st_held_control` | ⛔ **VOID** — `MNO_MOM=0` disabled the MOM that `MNO_ANNEAL_PENALTY` needs; both then SEGFAULTED (rc=139), flag-independent | — | — | — |

**Read Eee.**  It orders the failures by charge slosh (13.5 → 29.0 → 35.1) and separates the two distinct
mechanisms: the flag and flat-Kerker SLOSH, while linear-D-mix loses the moment with Eee normal.

## 2026-08-25 (cont.) — N2: is the ρ<0 ringing ALIASING or BAND-LIMITING?

| log | arm | pts ρ<0 (matched, end of stage 1) | neg mass (e) | min ρ | Etot | CPU / RSS |
|---|---|---|---|---|---|---|
| `st_baseline` | C=2, **BallOnly** (production) | 15.92% | −1.81e-3 | −0.146 | −61.40297621 | 500 s / 1316 MB |
| `st_aliasfree` | C=2, **AliasFree** — *same {G}*, 8× raster pts | **15.83%** | **−1.82e-3** | **−0.147** | −61.40297359 | 1155 s / 4652 MB |
| `st_C3` | C=3 (ecut 108) | 14.19% | −1.88e-4 | −0.0366 | −61.40282017 | 619 s / 1973 MB |
| `st_C4` | C=4 (ecut 144) | 6.55% | −1.65e-5 | −0.0071 | −61.40283283 | 814 s / 2638 MB |
| `st_C6` | C=6 (ecut 216), stopped after stage 1 | 5.61% | −3.02e-7 | −5.49e-4 | — | — |

**The raster A/B holds {G} FIXED** (nG 8623, 870 G-stars in both arms; BallOnly level-0 40³, AliasFree 80³).
8× the real-space points changes the lobes by **1%** ⇒ they are **Gibbs ringing from band-limiting, not
fold-back aliasing**.  Widening the ball (C) *does* kill them — ~6000× in negative mass by C=6 — because it
enlarges {G}, not because it reduces aliasing.

**Two verdicts fall out:** `BallOnly` is vindicated as the default (the never-taken A/B the code asks for —
fold-back is worth 2.6 µHa, AliasFree costs 2.3× CPU / 3.5× RAM); and since raster geometry does nothing
while a bigger ball does not scale (8× volume for a 2×2×2 supercell ⇒ >10 GB of raster at today's C=2),
the only remaining cure is the LOCAL cusp deficit — `doc/OpenWork.md` N4.
