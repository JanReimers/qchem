import sys, re, numpy as np
from pyscf.pbc import gto
BS='/home/janr/Code/qchem6/BasisSetData/valence_lowq_va.bsd'
def parse_bsd(fn):
    out={}; el=None
    for ln in open(fn):
        t=ln.split()
        if len(t)==2 and t[1]=='0' and t[0].isalpha(): el=t[0].capitalize(); out[el]=[]; continue
        if len(t)==3 and t[0] in 'SPDF' and el: l='SPDF'.index(t[0]); out[el].append([l,None]); continue
        if len(t)==2 and el and out[el] and out[el][-1][1] is None: out[el][-1][1]=float(t[0].replace('D','E'))
    return {e:[[l,[x,1.0]] for l,x in v] for e,v in out.items()}
def load(*fns):
    d={}
    for f in fns: d.update(parse_bsd(f))
    return d
A=0.529177210903
def rutile():
    a,c,u=4.594/A,2.959/A,0.305
    return [[a,0,0],[0,a,0],[0,0,c]],[('Ti',(0,0,0)),('Ti',(.5,.5,.5)),('O',(u,u,0)),('O',(-u,-u,0)),('O',(.5+u,.5-u,.5)),('O',(.5-u,.5+u,.5))]
def perov(a,M,A_,X):
    a=a/A
    return [[a,0,0],[0,a,0],[0,0,a]],[(A_,(0,0,0)),(M,(.5,.5,.5)),(X,(.5,.5,0)),(X,(.5,0,.5)),(X,(0,.5,.5))]
def rock(a,M,X,afm=True):
    a=a
    return [[0,a/2,a/2],[a/2,0,a/2],[a/2,a/2,0]],[(M,(0,0,0)),(X,(.5,.5,.5))]
def rockAFM(a,M,X):
    L=[[a,a/2,a/2],[a/2,a,a/2],[a/2,a/2,a]]  # 4 atom cell as in qchem (lattice a/2*(...) T-rebased)
    return L,[(M,(0,0,0)),(M,(.5,.5,.5)),(X,(.25,.25,.25)),(X,(.75,.75,.75))]
def pyrite(M='Fe',X='S'):
    a=5.417/A; u=0.385
    pos=[(0,0,0),(.5,.5,0),(.5,0,.5),(0,.5,.5)]
    S=[]
    for (x,y,z) in [(u,u,u),(-u+.5,-u,u+.5),(-u,u+.5,-u+.5),(u+.5,-u+.5,-u)]:
        S+= [(x,y,z),(-x,-y,-z)]
    return [[a,0,0],[0,a,0],[0,0,a]],[(M,p) for p in pos]+[(X,p) for p in S]
def k2np(a,c,z=0.153,zs=0.355):
    a=a/A; c=c/A
    return [[a,0,0],[0,a,0],[0,0,c]],[('Fe',(0,0,0)),('Fe',(.5,.5,.5)),('O',(.5,0,0)),('O',(0,.5,0)),('O',(0,.5,.5)),('O',(.5,0,.5)),
      ('O',(0,0,z)),('O',(0,0,-z)),('O',(.5,.5,.5+z)),('O',(.5,.5,.5-z)),('Sr',(0,0,zs)),('Sr',(0,0,-zs)),('Sr',(.5,.5,.5+zs)),('Sr',(.5,.5,.5-zs))]
def layered(M,a,c,zO=0.26):  # R-3m LiMO2 hex cell
    a=a/A;c=c/A
    L=[[a,0,0],[-a/2,a*np.sqrt(3)/2,0],[0,0,c]]
    atoms=[]
    for sh in [(0,0,0),(2/3,1/3,1/3),(1/3,2/3,2/3)]:
        atoms+= [(M,sh),('Li',tuple(np.add(sh,(0,0,.5)))),('O',tuple(np.add(sh,(0,0,zO)))),('O',tuple(np.add(sh,(0,0,-zO))))]
    return L,atoms
def S_min(basis,geom,kfr=((0,0,0),(.5,.5,.5),(.5,0,0),(.5,.5,0))):
    L,atoms=geom
    cell=gto.Cell(); cell.unit='B'; cell.a=np.array(L)
    cell.atom=[(el,tuple(np.dot(fr,L))) for el,fr in atoms]
    els=set(e for e,_ in atoms); cell.basis={e:basis[e] for e in els}
    cell.precision=1e-9; cell.verbose=0; cell.build()
    ks=cell.get_abs_kpts(np.array(kfr)); S=cell.pbc_intor('int1e_ovlp',kpts=ks)
    return cell.nao, min(np.linalg.eigvalsh(s)[0] for s in S)

def mindist(geom):
    L,at=geom; L=np.array(L); out={}
    import itertools
    for i,(e1,f1) in enumerate(at):
        for j,(e2,f2) in enumerate(at):
            for n in itertools.product(range(-1,2),repeat=3):
                d=np.linalg.norm((np.array(f2)+n-np.array(f1))@L)
                if d<1e-6: continue
                k=tuple(sorted((e1,e2))); out[k]=min(out.get(k,9e9),d)
    return {k:round(v,2) for k,v in out.items()}
