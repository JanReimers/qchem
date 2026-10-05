from cond import *
def G(M):
    d={'Ti':{'TiO2':rutile()},
     'V':{'SrVO3':perov(3.84,'V','Sr','O'),'LiVO2':layered('V',2.84,14.8)},
     'Cr':{'LiCrO2':layered('Cr',2.90,14.4)},
     'Fe':{'Sr2FeO4':k2np(3.86,12.4),'LiFeO2':layered('Fe',3.02,14.5),'FeS2':pyrite()},
     'Co':{'LiCoO2':layered('Co',2.816,14.05)},
     'Cu':{'KCuF3':perov(4.14,'Cu','K','F')},
     'Ni':{'NiO':rockAFM(7.88,'Ni','O'),'LiNiO2':layered('Ni',2.88,14.19)},
     'S':{'FeS2':pyrite()},'K':{'KCuF3':perov(4.14,'Cu','K','F')},'Sr':{'SrVO3':perov(3.84,'V','Sr','O'),'Sr2FeO4':k2np(3.86,12.4)}}
    return d[M]
if __name__=='__main__':
    for m in ['Ti','V','Fe','Co','Cu','Ni']:
        for k,g in G(m).items(): print(k,mindist(g))
