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
