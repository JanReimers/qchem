import numpy as np
from pyscf import gto
HA=27.211386245988
def F0(l,exps,cs,omega,elem):
    mol=gto.M(atom=f'{elem} 0 0 0',verbose=0,basis={elem:[[l]+[[e,c] for e,c in zip(exps,cs)]]})
    n=mol.nao
    with mol.with_range_coulomb(omega): eri=mol.intor('int2e').reshape(n,n,n,n)
    return np.mean([eri[m,m,p,p] for m in range(n) for p in range(n)])
NI=(2,[0.180,0.39485815,0.86618308,1.90010803,4.16818406,9.14356349,20.05783621,44.0],
      [-0.114,-0.154,-0.292,-0.308,-0.387,-0.047,0.009,-0.002],'Ni')
O =(1,[0.465,1.200,3.098,8.000],[0.705,0.092,0.335,0.045],'O')
def need(man,target):
    l,e,c,el=man; bare=F0(l,e,c,0.0,el)*HA; want=target/bare
    lo,hi=0.01,5.0
    for _ in range(60):
        m=0.5*(lo+hi)
        if F0(l,e,c,-m,el)*HA/bare>want: lo=m
        else: hi=m
    return 0.5*(lo+hi),bare
wNi,bNi=need(NI,5.27)
print(f"Ni 3d: bare {bNi:.2f} eV, target 5.27 (hp.x) -> omega {wNi:.3f}")
print("HISTORICAL SENSITIVITY SCAN -- O 2p used to be a BOUND (>~4); A2b (2026-09-25) replaced it with a")
print("real matched-PP hp.x value, 8.5139 eV (IntegrationTests/QE/README.md NiOgO.*), added to the scan below.")
for tgt in (3.0,4.0,5.0,6.0,7.0,7.27,8.5139):
    w,b=need(O,tgt)
    print(f"  O 2p target {tgt:4.2f} eV -> omega {w:.3f}   (disagreement with Ni 3d: {abs(w-wNi)/wNi*100:5.1f} %)")
print("\nAnd what a FLAT dielectric would need instead of a range cut:")
print(f"  Ni 3d: eps = {bNi/5.27:.2f}   (NiO experimental eps_inf ~ 5.7)")
