# Hamiltonian cleanup sprint

PLAN, born 2026-11-06 (user: geometry optimization is Hamiltonian centred, so clear the H rows of `doc/OOD-SOLID-Cleanup.md` first).  Under
execution after the Saturday 2026-11-07 reboot.  One lane in sequence, lane S: these rows share hot files.  Retires to `doc/Records/` when
step 5 lands.  The rows stay the record of what is open; this file holds only the order, the gates and the rulings not yet on a row.

## Order

| step | rows | why here | gate |
|---|---|---|---|
| 1 | `§I.1`, `D13`, `BM(1)` | independent, small | sweep; a move is re-judged against the oracle |
| 2 | `V1.34` | ruled shape (3) not yet realised, and it blocks the rest | sweep, 0 moves |
| 3 | `R1.0q`: `FittedVxc`, then `Vee_Hartree` and the quadrature XC terms | needs step 2 | sweep, 0 moves |
| 3b | `V1.1b`: `Vee_DFT`, Dunlap and normal | see the ruling below; touches the term step 3 passes through | design note first; sweep, 0 moves on the default route |
| 4 | `V1.35` with `V1.41` stages 2-3, planned together | the `*_RealBlock` caches are in both; both change the tier-1 `Hamiltonian` hub | design notes first; sweep + panel at 4 threads; one timed MnO or TeF₆ run before and after |
| 5 | `V-XC-ONE-TERM`, then `V-NONLINEAR-MESH` | the term loses mesh ownership in the second | `VxcProductSpace`; moves re-judged against PySCF; one timed run |

Not in the sprint: `V1.39` (its own trigger is GGA), `V-BLOCKED-MATRIX` (it may fall out of step 4; decide there).

**While in each file (user 2026-11-06, `OOD-SOLID-Cleanup.md` header):** concrete classes not seen by clients are low priority, but look for the same sequence of code appearing in two structure-specific classes, and log or fold it.

## How a step is validated

No row needs a material campaign.  Those came from numerical questions; a cleanup asks whether the numbers stayed the same.

- **0 moves expected:** short + long, `scripts/anchordiff --tol 1e-10`.  About 1.3 hours.
- **Inside the Fock loop** (`V1.41`, the `R1.0r` magnetic branch): also the panel at 4 threads.  Hours, unattended.
- **Moves on purpose** (`V-XC-ONE-TERM`): re-judged against PySCF on small molecules (pin 10).
- **Not seen by a sweep:** wall time, so one timed run where caches move; threads (baselines are serial); opt-in routes, which need an
  energy gate, not only a pointwise one (the radial Poisson bug, `b0f1ea7e`).

## Design notes the user sees first

`Hamiltonian` is a tier-1 hub, so these gate the sprint more than validation does.

1. **`V1.35`: the blast radius (user 2026-11-06).**  `tDynamic_HT<TBlock,TRun>` looks clean by itself.  The real question: do the two
   parameters propagate through existing code, or do many classes become templates that were not?  That may be fine, but see it first.
   The note counts it: which classes change, which become templates, which stay concrete because the parameters stop at a boundary, how
   many module interface units are touched, and what that does to build time.  `*RealBlock` is named in about 20 files today.
2. **`V1.41` stages 2-3:** the `MakeMatrix(bs, s, DensityID)` term interface (ruled), and what moves per term.
3. **`V1.1b`:** below.

## Ruling: `V1.1b` (user 2026-11-06, to be copied to its row)

The Dunlap fit is not named explicitly today, and it should be.  Target: two versions of the Hartree term (`FittedVee`, or `Vee_DFT` under
the new naming), one Dunlap and one normal.  Each creates the correct fitter implementations and calls the correct 3-centre integrals.
Today one `FittedVee` serves both and tells them apart by `itsGridCharge`, and the energy formula 2·EeeFit − EeeFitFit is valid only for
Dunlap: that is the defect.  The two versions are the shapes in `doc/Records/HartreeFitConsistency.md`: Dunlap Cᵀ G⁻¹ C, second order;
normal Oᵀ K O, first order, with `K` a `CoulombOperator`.

**Open for the note: can they be structure neutral?**  Wanted, if achievable.  The likely obstacle is G = 0: a periodic Coulomb metric
diverges for a charged fit function, so a periodic Dunlap fit needs its neutrality stated, where a normal fit in G-space drops G = 0 in `K`.
Say what each version does on a molecule, an atom and a solid, and where neutrality truly fails, name it.
