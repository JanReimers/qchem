// File: ElectronConfiguration.C

#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/Misc/PeriodicTable.H"
#include "Imp/Symmetry/AngularQN.H"
#include <Spin.H>
#include <cassert>
#include <iostream>
#include <initializer_list>

using std::cout;
using std::endl;

PeriodicTable pt;

const int Atom_EC::FullShells[Nshell][LMax+2]=
{ // Z   s p d f
    {0  ,0,0,0,0}, //1st row
    {2  ,1,0,0,0}, //2nd row
    {10 ,2,1,0,0},
    {18 ,3,2,0,0},
    {36 ,4,3,1,0},
    {54 ,5,4,2,0},
    {86 ,6,5,3,1},
    {118,7,6,4,2}
};

Atom_EC::Atom_EC(int Z)
: N{0,0,0,0}
, NUnpaired(pt.GetNumUnpairedElectrons(Z))
{
    assert(Z>0);
    assert(Z<=N_Elements);
    int* vc=pt.GetValanceConfiguration(Z);  //Valance electron counts
    // Hunt for the nearest for full shell for Z.
    int ns=0;
    for (;ns<Nshell;ns++)
        if (FullShells[ns][0]>Z) break;
    ns--;
//    for (auto l:{0,1,2,3}) cout << valance_configuration[l] << ",";
//    cout << endl;
//    for (auto l:{0,1,2,3}) cout << FullShells[ns][l+1] << ",";
//    cout << endl;

    // Load up arrays 
    //   Nf=# of {s,p,d,f} electrons in the full shells.
    //   Nv=# of {s,p,d,f} valance electrons
    for (int l=0;l<=LMax;l++) 
    {
        int g=2*(2*l+1); //degeneracy
        Nf[l]=FullShells[ns][l+1]*g;
        Nv[l]=vc[l];
        if (Nv[l]>0 && Nv[l]%g==0) //Check of valance shell is full.
        {   
            Nv[l]-=g;
            Nf[l]+=g;
        }
        N[l]=Nf[l]+Nv[l]; //Total electrons for s,p,d,f
        //cout << N[l] << ",";
    }
    //cout << endl;
}

int Atom_EC::GetN() const
{
    int ne=0;
    for (auto n:N) ne+=n; //Sum over l
    return ne;
}

int Atom_EC::GetN(const Spin& s) const
{
    if (s==Spin::None) return GetN();
    int ne=GetN();
    assert((ne+NUnpaired)%2==0);
    return s==Spin::Up ? (ne+NUnpaired)/2 : (ne-NUnpaired)/2;
}

int Atom_EC::GetN(const Symmetry& qn) const
{
    const AngularQN& sqn=dynamic_cast<const AngularQN&>(qn);
    int nl,nlu;
    std::tie(nl,nlu)=sqn.GetN(N,Nv,NUnpaired);
    return nl;    
}
int Atom_EC::GetN(const Symmetry& qn, const Spin& s) const
{
    if (s==Spin::None) return GetN(qn);
    
    const AngularQN& sqn=dynamic_cast<const AngularQN&>(qn);
    int nl,nlu;
    std::tie(nl,nlu)=sqn.GetN(N,Nv,NUnpaired);
    assert((nl+nlu)%2==0);
    return s==Spin::Up ? (nl+nlu)/2 : (nl-nlu)/2;            
}

void Atom_EC::Display() const
{
    cout << "N: ";
    for (auto n:N) cout << n << ",";
    cout << endl;
    cout << "Nf: ";
    for (auto n:Nf) cout << n << ",";
    cout << endl;
    cout << "Nv: ";
    for (auto n:Nv) cout << n << ",";
    cout << endl;
    cout << "NUnpaired: " << NUnpaired << endl;
}
    

int Molecule_EC::GetN(const Spin& s) const
{
    if (s==Spin::None) return GetN();
    if (Ne%2==0)
        return Ne/2;
    else
        return s==Spin::Up ? (Ne+1)/2 : (Ne-1)/2;
}

void Molecule_EC::Display() const
{
    cout << "Ne: " << Ne << endl;
}
