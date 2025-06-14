// File: Atom_EC.C

#include "Symmetry/Atom_EC.H"
#include "Common/PeriodicTable.H"
#include "Symmetry/Angular.H"
#include <Symmetry/Irrep_QNs.H>
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
    assert(Z>0);
    assert(Z<=N_Elements);

    int* vc=pt.GetValanceConfiguration(Z);  //Valance electron counts
    
    // Step 1: Hunt for the nearest for full shell for Z.
    int ns=0;
    for (;ns<Nshell;ns++)
        if (FullShells[ns][0]>Z) break;
    ns--;

    // Step 2: Load up counts of valance and core electrons.
    //   Nf=# of {s,p,d,f} electrons in the full shells.
    //   Nv=# of {s,p,d,f} valance electrons
    //   N=Nf+Nv
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
   
    size_t lv=0,l=0;
    for (auto nv:itsNs.Nv)  //Find highest l for valance electrons.
    {
        if (nv>0) lv=l;
        l++;
    }
    // Step 3:  Figure out where all the unpaired electrons land.  This is complicated by atoms like Cr
    //          with 6 unpaired electrons.  5 are d, and 1 is s.

    int nup=pt.GetNumUnpairedElectrons(Z); //Total # of unpaired electrons.
    int gv=2*lv+1; //Degeneracy of the valance shell.
    if (itsNs.Nv[lv]<=gv)
    { // <= half filled
        if (nup<=itsNs.Nv[lv])
        {
            itsNs.Nu[lv]=nup; //take all unapired in the valance shell.
            nup=0;
        }
        else
        {
            itsNs.Nu[lv]=itsNs.Nv[lv];  //Take as many as we can in the valance shell.
            nup-=itsNs.Nv[lv]; //Keep track of how many unpaired electrons are left over.
        }

    }   
    else
    { //> half filled
        int Nvu=2*gv-itsNs.Nv[lv]; //# of unpaird orbitals.
        assert(Nvu<=gv);
        if (nup<=Nvu)
        {
            itsNs.Nu[lv]=nup; //take all unapired in the valance shell.
            nup=0;
        }
        else
        {
            itsNs.Nu[lv]=Nvu;  //Take as many as we can in the valance shell.
            nup-=Nvu; //Keep track of how many unpaired electrons are left over.
        }

    }
    //  Put any remaining unpaired electrons into lower l
    if (nup>0)
    {
        assert(nup==1); //There are no cases where lower l gets more than one unpaired.
        assert(lv>0);
        for (int l=lv-1;l>=0 && nup>0;l--) //Go backwards throuhg ls
            if (itsNs.Nv[l]>0) //Did we find some valance electrons
            {
                assert(itsNs.Nv[l]==1); 
                itsNs.Nu[l]=nup;
                nup--;
            }
    }
    // Display();
    assert(nup==0); //By now all unpaired electrons should have gobbled up.        
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
    int NUnpaired=0;
    for (auto& nu:itsNs.Nu) NUnpaired+=nu;
    assert((ne+NUnpaired)%2==0);
    return s==Spin::Up ? (ne+NUnpaired)/2 : (ne-NUnpaired)/2;
}

int Atom_EC::GetN(const sym_t& qn) const
{
    const Angular_Sym* sqn=dynamic_cast<const Angular_Sym*>(qn.get());
    ElCounts_l ecl=sqn->GetN(itsNs);
    return ecl.N; // Should be total core+valance 
}
int Atom_EC::GetN(const Irrep_QNs& qns) const
{
    if (qns.ms==Spin::None) return GetN(qns.sym);
    
    const Angular_Sym* sqn=dynamic_cast<const Angular_Sym*>(qns.sym.get());
    ElCounts_l ecl=sqn->GetN(itsNs);
    assert((ecl.N+ecl.Nu)%2==0);
    return ecl.GetN(qns.ms);  // Should be total core+valance      
}

ml_Breakdown Atom_EC::GetBreadown(size_t l) const
{
    itsNs.DebugCheck(); //Check self consistency
    ml_Breakdown mls;
    size_t g=(2*l+1); //degenracy
    size_t Nunp=itsNs.Nu[l];
    size_t Npairs= (itsNs.Nv[l]-itsNs.Nu[l])/2; // For full shell systems this ends up being zero.
    size_t Nempty=g-Npairs-Nunp;
    if (Nempty==g) //Fix up for full shell systems.
    {
        Npairs=g;
        Nempty=0;
    }
    int ml=-(int)l;
    for (size_t i=0;i<Npairs;i++) mls.ml_paired    .push_back(ml++);
    for (size_t i=0;i<Nunp  ;i++) mls.ml_unpaired  .push_back(ml++);
    for (size_t i=0;i<Nempty;i++) mls.ml_unoccupied.push_back(ml++);

    // cout << "ml_paired    =" << mls.ml_paired << endl;
    // cout << "ml_unpaired  =" << mls.ml_unpaired << endl;
    // cout << "ml_unoccupied=" << mls.ml_unoccupied << endl;
    assert(ml==(int)(l+1));
    return mls;
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
    cout << "Nu: ";
    for (auto n:itsNs.Nu) cout << n << ",";
    cout << endl;
}
   