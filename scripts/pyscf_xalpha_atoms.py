#!/usr/bin/env python3
# The Slater-Xalpha atom ORACLE (OpenWork "A_PG_DFT asserts the SIGNED relative error", lane A 2026-10-31).
#
# The Xalpha atom gates (IntegrationTests/A_DFT.C: A_DFT_atom SG/SL, A_PG_DFT) used to compare against NIST's LSDA
# (Dirac + VWN) -- a DIFFERENT functional -- and passed only because the assertion was signed.  This script computes the
# independent reference (doc/Pins.md pin 10): PySCF with the SAME functional, E_x scaled by alpha/(2/3) (xc = "1.5*alpha*SLATER",
# the convention our Xalpha matches: B at alpha 0.78744 agrees to the basis gap), at the SAME per-Z Schwarz alpha
# (src/Common/Imp/PeriodicTable.C theSlaterAlpha; 0.70 past Ca), in a near-complete uncontracted basis, non-relativistic.
# Two states per open-shell atom, because they are DIFFERENT states (doc/SCFRobustnessPlan.md §C1-DEG):
#   integer = UKS aufbau (a pure M_L component -- what a molecular-path atom with a broken-symmetry fill converges to),
#   frac    = UKS with the degenerate open shell fractionally occupied (the spherical M_L ensemble -- the radial atom code's state).
# Run with the PySCF venv:   ~/Code/pyscf-env/bin/python scripts/pyscf_xalpha_atoms.py
import sys
from pyscf import gto, dft, scf

ALPHA = {1:0.77679, 2:0.77224, 3:0.79118, 4:0.79526, 5:0.78744, 6:0.77657, 7:0.76654, 8:0.76454, 9:0.75954, 10:0.73100,
         11:0.75110, 12:0.74942, 13:0.74797, 14:0.74521, 15:0.74309, 16:0.74270, 17:0.74183, 18:0.722, 19:0.7227, 20:0.7220}
SYM = {1:'H',2:'He',3:'Li',4:'Be',5:'B',7:'N',10:'Ne',11:'Na',18:'Ar',36:'Kr',37:'Rb',53:'I',54:'Xe'}
SPIN = {1:1,2:0,3:1,4:0,5:1,7:3,10:0,11:1,18:0,36:0,37:1,53:1,54:0}   # 2S of the polarized ground state

def run(Z, polarized, frac, basis):
    a = ALPHA.get(Z, 0.70)
    m = gto.M(atom=f'{SYM[Z]} 0 0 0', basis=basis, spin=SPIN[Z] if polarized else 0, verbose=0)
    mf = (dft.UKS(m) if polarized else dft.RKS(m))
    mf.xc = f'{1.5*a}*SLATER'; mf.grids.level = 7; mf.conv_tol = 1e-11
    if frac: mf = scf.addons.frac_occ(mf)
    E = mf.kernel()
    if not mf.converged: E = float('nan')
    return a, E

cases = [(Z,False) for Z in (2,4,10,18,36,54)] + [(Z,True) for Z in (1,3,5,7,11,37,53)]
basis = sys.argv[1] if len(sys.argv)>1 else 'unc-ano-rcc'
for Z,pol in cases:
    a, Ei = run(Z, pol, False, basis)
    line = f"Z={Z:2d} {SYM[Z]:2s} {'pol  ' if pol else 'unpol'} alpha={a:.5f} {basis} integer E={Ei:.8f}"
    if pol and Z in (5,53):
        _, Ef = run(Z, pol, True, basis)
        line += f"  frac(spherical) E={Ef:.8f}"
    print(line, flush=True)
