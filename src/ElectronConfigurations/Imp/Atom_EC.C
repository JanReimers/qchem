// File: ElectronConfigurations/Imp/Atom_EC.C  Electron configuration for atoms.
module;
#include <cassert>
#include <iostream>
#include <initializer_list>
#include <vector>

module qchem.ElectronConfiguration.AtomNR;
import qchem.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Spherical;
import qchem.Symmetry.Factory;
import qchem.stl_io;
import qchem.Blaze;

using std::cout;
using std::endl;


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

Atom_EC::Atom_EC(int Z, NsOnly_t)
: itsLMax(0), itsLValance(0)
{
    assert(Z>0);
    assert(Z<=N_Elements);

    const size_t* vc=thePeriodicTable().GetValanceConfiguration(Z);  //Valance electron counts

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
    // Step 3:  Figure out where all the unpaired electrons land.
    AssignUnpaired(Z, thePeriodicTable().GetNumUnpairedElectrons(Z));
}

// Distribute the atom's unpaired electrons into itsNs.Nu.  This is complicated by atoms like Cr with 6
// unpaired electrons (5 are d, 1 is s), so a single leftover spills into the next-lower partially-filled l.
void Atom_EC::AssignUnpaired(int Z, int nup)
{
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
    assert(nup==0); //By now all unpaired electrons should have been gobbled up.
}

void Atom_EC::SetSplitOccupations(sym_t sp, sym_t su, int NCore, int gp, int gu, int Npair, int Nu)
{
    itsOccupations[Irrep(Spin::Up  ,sp)]=NCore*gp+Npair;
    itsOccupations[Irrep(Spin::Down,sp)]=NCore*gp+Npair;
    itsOccupations[Irrep(Spin::Up  ,su)]=NCore*gu+Nu;
    itsOccupations[Irrep(Spin::Down,su)]=NCore*gu;
}

void Atom_EC::BuildNROccupations()
{
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
            // Npair is the partially-filled-but-PAIRED part (nonzero only for an un-demoted closed
            // valence subshell, e.g. the pseudo-atom's Si 3s^2); for normal atoms it is 0 here because
            // closed valence subshells fold into Nf (NCore).  Nu is nonzero only in the gp==0 half-filled case.
            sym_t s=Symmetry::YFactory(l);
            itsOccupations[Irrep(Spin::Up  ,s)]=NCore*g+Npair+Nu;
            itsOccupations[Irrep(Spin::Down,s)]=NCore*g+Npair;
        }
        else //ml splitting: unpaired electrons occupy the lowest ml values
        {
            assert(Nv<2*g);
            ivec_t ms_u(gu), ms_p(gp);
            int ml=-(int)l;
            for (int& m:ms_u) m=ml++;
            for (int& m:ms_p) m=ml++;
            sym_t su=Symmetry::YFactory(l,ms_u);
            sym_t sp=Symmetry::YFactory(l,ms_p);
            SetSplitOccupations(sp,su,NCore,gp,gu,Npair,Nu);
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
}

Atom_EC::Atom_EC(int Z)
: Atom_EC(Z, NsOnly_t{})
{
    BuildNROccupations();
}

// Pseudo-atom: load ONLY the raw valence shells (no core, no full-subshell demotion), so a closed valence
// subshell like Si 3s^2 stays a valence occupation instead of folding into the core.  netCharge != 0 makes
// an ion (see the header): adjust the neutral valence by -netCharge electrons, then recompute the unpaired
// count by Hund's rule on the (adjusted) valence shell.
Atom_EC::Atom_EC(int Z, int netCharge, ValenceOnly_t)
: itsLMax(0), itsLValance(0)
{
    assert(Z>0);
    assert(Z<=N_Elements);

    const size_t* vc=thePeriodicTable().GetValanceConfiguration(Z);  //Raw NEUTRAL valence electron counts per l
    int nv[LMax+1];
    for (size_t l=0;l<=LMax;l++) nv[l]=int(vc[l]);

    // Make the ion: add (anion) / remove (cation) -netCharge electrons.  Anion fills the lowest l that has
    // room (capacity 2(2l+1)); cation removes from the highest l that still has electrons.  Main-group ions
    // (the IonicSAD targets) close or empty a shell this way -- F- -> p^6, Na+ -> empty, O2- -> p^6.
    for (int add = -netCharge; add>0; --add)         // anion: add electrons low-l first
    {
        int l=0; while (l<=int(LMax) && nv[l]>=2*(2*l+1)) ++l;
        assert(l<=int(LMax) && "PseudoAtom_EC: anion overflows the valence (no open channel <= LMax)");
        ++nv[l];
    }
    for (int rem = netCharge; rem>0; --rem)          // cation: remove electrons high-l first
    {
        int l=int(LMax); while (l>=0 && nv[l]<=0) --l;
        assert(l>=0 && "PseudoAtom_EC: cation removes more electrons than the valence holds");
        --nv[l];
    }

    for (size_t l=0;l<=LMax;l++)
    {
        itsNs.Nf[l]=0;                 //no core: the pseudopotential replaces it
        itsNs.Nv[l]=nv[l];
        itsNs.N [l]=nv[l];
        if (nv[l]>0) {itsLMax=l; itsLValance=l;}
    }
    // Neutral atom: keep the periodic table's exact unpaired count (handles half-filled-shell exceptions like
    // Cr/Cu).  Ion: recompute by Hund's rule on the (highest-l) adjusted valence shell -- n in 2g spin-orbitals
    // gives n unpaired if n<=g else 2g-n, so a closed-shell ion (F-/Na+/O2-) is 0.  (Open-shell transition-
    // metal ions are out of scope.)
    int nup;
    if (netCharge==0)
        nup = (int)thePeriodicTable().GetNumUnpairedElectrons(Z);
    else
    {
        int g=2*int(itsLValance)+1, nlv=nv[itsLValance];
        nup = (nlv<=g) ? nlv : (2*g-nlv);
    }
    AssignUnpaired(Z, nup);
}

PseudoAtom_EC::PseudoAtom_EC(int Z, int netCharge)
: Atom_EC(Z, netCharge, ValenceOnly_t{})
{
    BuildNROccupations();
}

int Atom_EC::GetN(const Irrep& qns) const
{
    if (qns.ms==Spin::None)
    {
        // A basis irrep that has no occupation entry is a legitimately-EMPTY (virtual) irrep -- e.g. the
        // m=+1 p sub-irrep of a pseudo-atom's open p^2 shell, which is occupied only in p{-1,0}.  Return 0.
        // (The end()-comparison was previously against the WRONG container, masking the miss with UB.)
        auto i=itsUnpolOccupations.find(qns);
        return i==itsUnpolOccupations.end() ? 0 : i->second;
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
   
