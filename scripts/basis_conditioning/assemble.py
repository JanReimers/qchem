from scan import *
W={'Li':(1,['0:4:0.08:2']),'Ti':(4,['0:7:0.10:24','2:7:0.18:36']),'V':(5,['0:7:0.10:24','2:7:0.25:36']),
 'Cr':(6,['0:7:0.10:24','2:7:0.25:36']),'Fe':(8,['0:7:0.10:24','2:7:0.25:36']),'Co':(9,['0:7:0.10:24','2:7:0.30:36']),
 'Ni':(10,['0:7:0.13:24','2:7:0.30:36']),'Cu':(11,['0:7:0.10:24','2:7:0.25:36']),'K':(1,['0:6:0.06:8']),
 'Sr':(2,['0:6:0.06:10']),'S':(6,['0:6:0.2:30','1:6:0.16:20'])}
blocks={}; gaps={}
for e,(q,sh) in W.items():
    b,g,fn=mk(e,q,sh,'final'); blocks[e]=b; gaps[e]=g
for e in ('O','F'): blocks[e]=base[e]
def L(l): return 'SPDF'[l]
bsd='''valence_lowq_tm: GTH-pseudopotential valence bases for the elements of the Hubbard-U oracle set, minted 2026-10-04
for the CP2K finite-alpha LRT comparison (doc/HubbardUPlan.md) and conditioned for the SOLID, not the free atom.
Each window is the most diffuse one (smallest valgen gap to the complete-basis PP floor) whose overlap S(k)
keeps lambda_min >= 4e-4 on every host crystal of the element (idealised geometries: rutile TiO2, cubic SrVO3 and
KCuF3, K2NiF4 Sr2FeO4, pyrite FeS2, rocksalt NiO, R-3m LiMO2; k = G, X, L, (1/2,1/2,0)); the accepted Mn VA basis
is 6.4e-4 on MnO.  Measured with scratch tooling (PySCF pbc overlap), NOT with the vet stage.
Lessons: (1) the free-atom gap is a poor judge -- NiO with the committed Ni (s .06, d .18) has lambda_min 1e-6;
(2) Li s 0.03 (the valence_lowq_va Li) alone drove LiMO2 to 1e-7, and TRIMMING a short even-tempered set
makes it WORSE (5 exps 0.5..2 -> 7e-6: beta ~1.4 is internally dependent) -- cut N and keep beta ~2.2;
(3) in the 3d series the diffuse s hurts conditioning more than the diffuse d.
Per-element recipe (valgen --q Q --shell ... --floor), gap = mHa above the pool floor:
''' + '\n'.join(f'  {e:2s} q{W[e][0]:<2d} {" ".join("--shell "+s for s in W[e][1]):45s} gap {gaps[e]:6.1f}' for e in W) + '''
O and F are the valence_lowq_va blocks (already diffuse-trimmed), repeated so this file is self-contained.

!
BASIS="Low-q GTH valence basis TM (generated, bulk-conditioned)"
'''
for e,b in blocks.items():
    bsd+=f' {e.upper()}   0\n'
    for l,x in [(a[0],a[1][0]) for a in b]: bsd+=f' {L(l)}   1  1.00\n            {x:.8f}        1.00000000\n'
    bsd+=' ****\n'
open('valence_lowq_tm.bsd','w').write(bsd)
cp=''
for e,b in blocks.items():
    cp+=f'{e} VALENCE-LOWQ-TM\n {len(b)}\n'
    for l,x in [(a[0],a[1][0]) for a in b]: cp+=f' {max(l+1,2)} {l} {l} 1 1\n {x:14.8f}  1.0\n'
open('VALENCE-LOWQ-TM-BASIS','w').write(cp)
print({e:len(b) for e,b in blocks.items()})
