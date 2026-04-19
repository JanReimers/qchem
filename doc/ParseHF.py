import os
import json
# print(os.getcwd()) where is the IDE running this?

#
# Script to parse the the data from Saito paper S.L. Saito / Atomic Data and Nuclear Data Tables 95 (2009) 836–870 
# Do: pdftotext -layout saito.pdf
# This script the parses saito.txt and output saito.json
#
class Orbital:
    def __init__(self,name,e,rs):
        self.name=name
        self.e=float(e)
        self.rs=rs # <r^2> <r> <1/r> <1/r^2> (for non-s orbitals <1/r^3> )
    def show(self):
        print("  ",self.name," e=",self.e)

class Element:
    def __init__(self,Z, symbol, valance,term):
        self.Z = Z
        self.symbol = symbol
        self.valance = valance
        self.term = term
        self.Orbitals=[]
    def SetHF(self,e):
        self.HFEnergy=float(e)
    def Setq0(self,q0):
        self.q0=float(q0)
    def AddOrbital(self,orbital):
        self.Orbitals.append(orbital)
    def show(self):
        print(self.symbol," Z=",self.Z," ",self.valance," ",self.term," E=",self.HFEnergy," q0=",self.q0)
        for o in self.Orbitals:
            o.show()
    

with open("doc/saito.txt", mode="r", encoding="utf-8") as file:
    inOrbitals=False
    el=None
    Elements=[]
    for line in file:
        if "Z = " in line:
            # found a new element.
            if el!=None:
                Elements.append(el);
            inOrbitals=False
            s1=line.split(',')
            symbol=s1[0]
            Z=int(s1[1][5:])
            s2=s1[2].split('(')
            valance=s2[0]
            term=s2[1][:-2]
            el=Element(Z,symbol,valance,term)
        elif line[:15]=="Total energy = ":
            print(line)
            HFenergy=line.split(';')[0].split('=')[1]
            el.SetHF(HFenergy)
        elif line[:5]=="q0 = ":
            q0=line.split(';')[0].split('=')[1]
            el.Setq0(q0)
        #     # print("q0=",q0)
        elif "EORB" in line:
            #found the begining of the orbitals list.
            inOrbitals=True
        elif inOrbitals==True:
            s1=line.split()
            if len(s1)==0 or s1[1]=="Saito" or s1[2]=="Saito": #Junk lines
                inOrbitals=False
            elif s1[0]=="(continued": #Junk lines
                inOrbitals=True
            else:
                # found and orbital record.
                el.AddOrbital(Orbital(s1[0],s1[1],s1[2:]))

            
for e in Elements:
    e.show()

with open("doc/saito.json", mode="w", encoding="utf-8") as write_file:
    json.dump(Elements, write_file, default=vars)