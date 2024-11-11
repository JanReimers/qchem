// File: ElectronConfiguration.C

#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/Misc/PeriodicTable.H"
#include <Spin.H>
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include <cassert>
#include <iostream>
#include <initializer_list>

using std::cout;
using std::endl;

PeriodicTable pt;

const int AtomElectronConfiguration::FullShells[Nshell][LMax+2]=
{
    {0  ,0,0,0,0},
    {2  ,1,0,0,0},
    {10 ,2,1,0,0},
    {18 ,3,2,0,0},
    {36 ,4,3,1,0},
    {54 ,5,4,2,0},
    {86 ,6,5,3,1},
    {118,7,6,4,2}
};

AtomElectronConfiguration::AtomElectronConfiguration(int Z)
: N{0,0,0,0}
, NUnpaired(pt.GetNumUnpairedElectrons(Z))
{
    assert(Z>0);
    assert(Z<=N_Elements);
    int* vc=pt.GetValanceConfiguration(Z);
    int ns=0;
    for (;ns<Nshell;ns++)
        if (FullShells[ns][0]>Z) break;
    ns--;
//    for (auto l:{0,1,2,3}) cout << valance_configuration[l] << ",";
//    cout << endl;
//    for (auto l:{0,1,2,3}) cout << FullShells[ns][l+1] << ",";
//    cout << endl;
    for (int l=0;l<=LMax;l++) 
    {
        Nf[l]=FullShells[ns][l+1]*2*(2*l+1);
        Nv[l]=vc[l];
        int g=2*(2*l+1); //degeneracy
        if (Nv[l]>0 && Nv[l]%g==0)
        {   
            Nv[l]-=g;
            Nf[l]+=g;
        }
        N[l]=Nf[l]+Nv[l];
        //cout << N[l] << ",";
    }
    //cout << endl;
}

int AtomElectronConfiguration::GetN() const
{
    int ne=0;
    for (auto n:N) ne+=n;
    return ne;
}

int AtomElectronConfiguration::GetN(const Spin& s) const
{
    if (s==Spin::None) return GetN();
    int ne=GetN();
    assert((ne+NUnpaired)%2==0);
    assert(s!=Spin::None);
    return s==Spin::Up ? (ne+NUnpaired)/2 : (ne-NUnpaired)/2;
}

int AtomElectronConfiguration::GetN(const QuantumNumber& qn) const
{
    const SphericalSymmetryQN& sqn=dynamic_cast<const SphericalSymmetryQN&>(qn);
    return N[sqn.GetL()];    
}
int AtomElectronConfiguration::GetN(const QuantumNumber& qn, const Spin& s) const
{
    if (s==Spin::None) return GetN(qn);
    const SphericalSymmetryQN& sqn=dynamic_cast<const SphericalSymmetryQN&>(qn);
    int l=sqn.GetL();
    int nl=N[l];
    if (Nv[l]==0) return nl/2;
    assert(nl!=0);
    // Handle partial shells
    int nlu=1; //# unpaired in shell l. 
    if (l==1) // p is partial.
    {
        assert(Nv[2]==0); //No partial D orbital
        nlu=NUnpaired-Nv[0];
    }
    else if (l==2) // d is partial.
    {
        assert(Nv[1]==0); //p better be full
        if (Nv[l]>1) nlu=NUnpaired-Nv[0];            
    }
    else if(l==3) // f is partial.
    {
        
        assert(Nv[0]==0); //If f is Partial s must be full.
        assert(Nv[1]==0); //If f is Partial p must be full.
        nlu=NUnpaired-Nv[2];
    }
    assert((nl+nlu)%2==0);
    return s==Spin::Up ? (nl+nlu)/2 : (nl-nlu)/2;            
}
    

int MoleculeElectronConfiguration::GetN(const Spin& s) const
{
    if (Ne%2==0)
        return Ne/2;
    else
        return s==Spin::Up ? (Ne+1)/2 : (Ne-1)/2;
}
