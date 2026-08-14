# Build probe .bsd variants from the committed v2 sph file by REMOVING added groups.
import subprocess, sys
v2 = subprocess.run(['git','-C','/home/janr/Code/qchem6','show','HEAD:BasisSetData/valence_lowq_sph.bsd'],
                    capture_output=True, text=True).stdout
def drop(t, block):
    assert block in t, block[:40]
    return t.replace(block, '', 1)
MND018 = ' D   1  1.00\n            0.18000000        1.00000000\n'
OS015  = ' S   1  1.00\n            0.15000000        1.00000000\n'
OP018  = ' P   1  1.00\n            0.18000000        1.00000000\n'
variant = sys.argv[1]
t = v2
if variant == 'A':            # SR + Mn7s ONLY (drop Mn d 0.18 + O diffuse)
    t = drop(t, MND018); t = drop(t, OS015); t = drop(t, OP018)
elif variant == 'B':          # SR + Mn7s + Mn d 0.18 (drop O diffuse)
    t = drop(t, OS015); t = drop(t, OP018)
elif variant == 'D':          # SR + O diffuse ONLY (drop Mn d 0.18; restore Mn 2s trim)
    t = drop(t, MND018)
    seven = ''' S   1  1.00
            0.10000000        1.00000000
 S   1  1.00
            0.24928829        1.00000000
 S   1  1.00
            0.62144650        1.00000000
 S   1  1.00
            1.54919334        1.00000000
 S   1  1.00
            3.86195754        1.00000000
 S   1  1.00
            9.62740780        1.00000000
 S   1  1.00
           24.00000000        1.00000000
'''
    two = ''' S   1  1.00
            0.10000000        1.00000000
 S   1  1.00
           24.00000000        1.00000000
'''
    assert seven in t
    t = t.replace(seven, two, 1)
open('/home/janr/Code/qchem6/BasisSetData/valence_lowq_sph.bsd','w').write(t)
print('variant', variant, 'written')
