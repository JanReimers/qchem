import subprocess, itertools, sys, os, re
from cond import *; from geoms import G
V='/home/janr/Code/qchem6/build/Release/CLIapps/valgen'
base=load(BS)
def mk(el,q,shells,tag):
    fn=f'{el}_{tag}.bsd'
    args=[V,'--element',el,'--q',str(q),'--floor','--nmax','100','--out',fn]
    for s in shells: args+=['--shell',s]
    o=subprocess.run(args,capture_output=True,text=True).stdout
    conv='converged   true' in o; m=re.search(r'gap_mHa\s+(\S+)',o)
    if not (conv and m) or not os.path.exists(fn): return None,None,None
    return parse_bsd(fn)[el.capitalize()],float(m.group(1)),fn
def tobasis(b): return [[l,[x,1.0]] for l,x in b]
