// File: Atom_ECImp.C  Electron configuration for atoms.
module;
#include <cassert>
#include <iostream>
#include <initializer_list>
#include <vector>

module qchem.Symmetry.AtomEC;
import Common.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Spherical;
import qchem.Symmetry.Factory;
import qchem.stl_io;

using std::cout;
using std::endl;

PeriodicTableSaito pt;

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
: itsLMax(0), itsLValance(0)
{
    assert(Z>0);
    assert(Z<=N_Elements);

    const size_t* vc=pt.GetValanceConfiguration(Z);  //Valance electron counts
    
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
   
    size_t l=0;
    for (auto nv:itsNs.Nv)  //Find highest l for valance electrons.
    {
        if (nv>0) itsLValance=l;
        l++;
    }
    // Step 3:  Figure out where all the unpaired electrons land.  This is complicated by atoms like Cr
    //          with 6 unpaired electrons.  5 are d, and 1 is s.

    int nup=pt.GetNumUnpairedElectrons(Z); //Total # of unpaired electrons.
    int gv=2*itsLValance+1; //Degeneracy of the valance shell.
    if (itsNs.Nv[itsLValance]<=gv)
    { // <= half filled
        if (nup<=itsNs.Nv[itsLValance])
        {
            itsNs.Nu[itsLValance]=nup; //take all unapired in the valance shell.
            nup=0;
        }
        else
        {
            itsNs.Nu[itsLValance]=itsNs.Nv[itsLValance];  //Take as many as we can in the valance shell.
            nup-=itsNs.Nv[itsLValance]; //Keep track of how many unpaired electrons are left over.
        }

    }   
    else
    { //> half filled
        int Nvu=2*gv-itsNs.Nv[itsLValance]; //# of unpaird orbitals.
        assert(Nvu<=gv);
        if (nup<=Nvu)
        {
            itsNs.Nu[itsLValance]=nup; //take all unapired in the valance shell.
            nup=0;
        }
        else
        {
            itsNs.Nu[itsLValance]=Nvu;  //Take as many as we can in the valance shell.
            nup-=Nvu; //Keep track of how many unpaired electrons are left over.
        }

    }
    //  Put any remaining unpaired electrons into lower l
    if (nup>0)
    {
        assert(nup==1); //There are no cases where lower l gets more than one unpaired.
        assert(itsLValance>0);
        for (int l=itsLValance-1;l>=0 && nup>0;l--) //Go backwards through ls
            if (itsNs.Nv[l]>0) //Did we find some valance electrons
            {
                assert(itsNs.Nv[l]==1); 
                itsNs.Nu[l]=nup;
                nup--;
            }
    }
    if (Z==58 && nup==0) //Special handling for Ce with 5d1↑ 4f1↓ valance configuration.
    {
        itsNs.Nu[2]=1; //One unpaired in d↑
        itsNs.Nu[3]=-1; //One unparied in f↓
    }
    //
    //  Now build a list of symmetries with occupation numbers.
    //
    for (int l=0;l<=LMax;l++)
    {
        int Nf=itsNs.Nf[l], Nv=itsNs.Nv[l], Nu=abs(itsNs.Nu[l]);
        int g=2*l+1; //degeneracy for one spin state
        assert(Nf%2==0);
        assert(Nf/2%g==0);
        assert(Nu<=g);
        assert((Nv-Nu)%2==0);
        int NCore=Nf/2/g;
        int Npair=(Nv-Nu)/2;
        int gu=Nu, gp=g-gu; //unpaired and paired degeneracies
        if (gu==0 || gp==0) //No ml splitting: shell is fully paired or half-filled
        {
            sym_t s=Symmetry::YFactory(l);
            itsOccupations[Irrep(Spin::Up  ,s)]=NCore*g+Nu;
            itsOccupations[Irrep(Spin::Down,s)]=NCore*g;
        }
        else //ml splitting: unpaired electrons occupy the lowest ml values
        {
            assert(Nv<2*g);
            ivec_t ms_u(gu), ms_p(gp);
            int ml=-(int)l;
            for (int i=0;i<gu;i++) ms_u[i]=ml++;
            for (int i=0;i<gp;i++) ms_p[i]=ml++;
            sym_t su=Symmetry::YFactory(l,ms_u);
            sym_t sp=Symmetry::YFactory(l,ms_p);
            itsOccupations[Irrep(Spin::Up  ,sp)]=NCore*gp+Npair;
            itsOccupations[Irrep(Spin::Down,sp)]=NCore*gp+Npair;
            itsOccupations[Irrep(Spin::Up  ,su)]=NCore*gu+Nu;
            itsOccupations[Irrep(Spin::Down,su)]=NCore*gu;
        }
    }
    //
    //  Now build an un polarized version of the itsOccupations list;
    //
    for (auto ir:GetIrreps())
    {
        Irrep nqns(Spin::None,ir), uqns(Spin::Up,ir),dqns(Spin::Down,ir);
        itsUnpolOccupations[nqns]=GetN(uqns)+GetN(dqns);
    }
        
   
    assert(nup==0); //By now all unpaired electrons should have gobbled up.        
}

int Atom_EC::GetN(const Irrep& qns) const
{
    if (qns.ms==Spin::None) 
    {
        auto i=itsUnpolOccupations.find(qns);
        if (i==itsOccupations.end())
        {
            std::cout << "Cannot find irrep=" <<  qns << endl;
            Display();
            exit(-1);
        }
        return i->second;
    }
    auto i=itsOccupations.find(qns);
    if (i==itsOccupations.end())
    {
        std::cout << "Cannot find irrep=" <<  qns << " seqn=" 
        << qns.SequenceIndex() << endl;
        Display();
        
        //Still need this for Dirac atoms.
        // const Spherical_Sym* sqn=dynamic_cast<const Spherical_Sym*>(qns.sym.get());
        // ElCounts_l ecl=sqn->GetN(itsNs);
        // assert((ecl.N+ecl.Nu)%2==0);
        // cout << qns << " N=" << ecl.GetN(qns.ms) << endl;
        return -1;  // Should be total core+valance    
    }
    // assert(i!=itsOccupations.end());
    return i->second;
     
}

// Get a list of spatial symmetries (ignore spin).  Use std::set to avoid duplicates.
Atom_EC::syms_t Atom_EC::GetIrreps() const
{
    syms_t syms;
    for (auto ir:itsOccupations)
        if (ir.second>0) 
            syms.insert(ir.first.sym);
    return syms;
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
    for (auto o:itsOccupations)
    {
        if (o.second!=0)
            cout << "N(" << o.first << ")=" << o.second << endl;
    }
    cout << "-------------------------------------------" << endl;
}
   
