#!/bin/bash
# TEMPLATE for a batch job: copy to ~/Code/qchem6-runs/batch/queue/NN_<name>.sh (the numeric prefix sets the order;
# run_queue.sh runs them one at a time).  The ONE rule (D-RUNDATA): the job writes its results where they will live,
# via scripts/rundir -- never into batch/work/<job> -- so nothing needs graduating afterwards.
#
#     ~/Code/<app>-runs/<Material>/<name>.<ext>        <app> = cp2k | qe | abinit | qchem6
#
# Prefix every output with what it is (mno_ckalpha_a0.out, not run1.out): several jobs share a Material directory.
# Put scratch (QE tmp/, wavefunction dumps) under tmp/ -- `scripts/rundir --audit` reports it by size.
set -euo pipefail
WORK=$(/home/janr/Code/qchem6/scripts/rundir qe FeS2)      # <- app, Material
cd "$WORK"

# ... generate decks, run the code (MPI codes ALWAYS via mpirun -- CLAUDE.md), e.g.
# mpirun -np 1 /home/janr/Code/q-e/PW/src/pw.x -i fes2.scf.in > fes2.scf.out
