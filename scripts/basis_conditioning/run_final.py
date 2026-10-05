from scan import *
THR=float(sys.argv[2]) if len(sys.argv)>2 else 3e-4
for e,q,sh in (('Li',1,['0:4:0.08:2']),('K',1,['0:6:0.06:8']),('Sr',2,['0:6:0.06:10']),('S',6,['0:6:0.2:30','1:6:0.16:20'])):
    base[e]=mk(e,q,sh,'chosen')[0]
el,q=sys.argv[1].split(':'); q=int(q); mats=G(el)
best=None
for se in (0.06,0.08,0.10,0.13):
  for de in (0.18,0.25,0.30,0.38,0.50):
    sh=[f'0:7:{se}:24',f'2:7:{de}:36']
    b,gap,fn=mk(el,q,sh,'x')
    if b is None: continue
    bas=dict(base); bas[el]=b
    mins={m:S_min(bas,g)[1] for m,g in mats.items()}
    ok=min(mins.values())>=THR
    if ok and (best is None or gap<best[0]): best=(gap,sh,mins)
print(el,q,'BEST' if best else 'NONE',best,flush=True)
