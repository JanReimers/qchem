// File: Atom_EC.C

#include "Imp/BasisSet/Atom/EC.H"
#include "Imp/Misc/PeriodicTable.H"
#include "Imp/BasisSet/Atom/Angular.H"
#include <Orbital_QNs.H>
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
: itsLMax(0)
{
    itsNs.NUnpaired=pt.GetNumUnpairedElectrons(Z);
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
    for (size_t l=0;l<=LMax;l++) 
    {
        int g=2*(2*l+1); //degeneracy
        itsNs.Nf[l]=FullShells[ns][l+1]*g;
        itsNs.Nv[l]=vc[l];
        if (itsNs.Nv[l]>0 && itsNs.Nv[l]%g==0) //Check of valance shell is full.
        {   
            itsNs.Nv[l]-=g;
            itsNs.Nf[l]+=g;
        }
        itsNs.N[l]=itsNs.Nf[l]+itsNs.Nv[l]; //Total electrons for s,p,d,f
        //cout << N[l] << ",";
        if (l>itsLMax && itsNs.N[l]>0) itsLMax=l;
    }
    //cout << endl;
}

int Atom_EC::GetN() const
{
    int ne=0;
    for (auto n:itsNs.N) ne+=n; //Sum over l
    return ne;
}

int Atom_EC::GetN(const Spin& s) const
{
    if (s==Spin::None) return GetN();
    int ne=GetN();
    assert((ne+itsNs.NUnpaired)%2==0);
    return s==Spin::Up ? (ne+itsNs.NUnpaired)/2 : (ne-itsNs.NUnpaired)/2;
}

int Atom_EC::GetN(const Symmetry& qn) const
{
    const Angular_Sym& sqn=dynamic_cast<const Angular_Sym&>(qn);
    ElCounts_l ecl=sqn.GetN(itsNs);
    return ecl.n;    
}
int Atom_EC::GetN(const Irrep_QNs& qns) const
{
    if (qns.ms==Spin::None) return GetN(*qns.sym);
    
    const Angular_Sym& sqn=dynamic_cast<const Angular_Sym&>(*qns.sym);
    ElCounts_l ecl=sqn.GetN(itsNs);
    assert((ecl.n+ecl.nu)%2==0);
    return ecl.GetN(qns.ms);       
}

void Atom_EC::Display() const
{
    cout << "N: ";
    for (auto n:itsNs.N) cout << n << ",";
    cout << endl;
    cout << "Nf: ";
    for (auto n:itsNs.Nf) cout << n << ",";
    cout << endl;
    cout << "Nv: ";
    for (auto n:itsNs.Nv) cout << n << ",";
    cout << endl;
    cout << "NUnpaired: " << itsNs.NUnpaired << endl;
}
   