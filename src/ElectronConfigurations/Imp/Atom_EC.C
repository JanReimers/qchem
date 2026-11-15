// File: ElectronConfigurations/Imp/Atom_EC.C  Electron configuration for atoms.
module;
#include <cassert>
#include <iostream>
#include <initializer_list>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <map>
#include <sstream>
#include <cctype>
#include <cstring>

module qchem.ElectronConfiguration.AtomNR;
import qchem.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Factory;
import qchem.Blaze;

namespace qchem {

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

//---------------------------------------------------------------------------------------------------
//  A STATED CONFIGURATION (2026-11-14, OpenWork row "All-electron ION seeds").
namespace {
const char kLetters[]="spdf";
int SubshellCapacity(int l) {return 2*(2*l+1);}
//! The noble-gas cores NIST's notation writes as [X], expanded.
const std::map<std::string, ShellOccupations>& NobleGasCores()
{
    static const std::map<std::string, ShellOccupations> cores=[]{
        std::map<std::string, ShellOccupations> c;
        ShellOccupations s={{1,0,2}};                                   c["He"]=s;
        s.insert(s.end(), {{2,0,2},{2,1,6}});                           c["Ne"]=s;
        s.insert(s.end(), {{3,0,2},{3,1,6}});                           c["Ar"]=s;
        s.insert(s.end(), {{3,2,10},{4,0,2},{4,1,6}});                  c["Kr"]=s;
        s.insert(s.end(), {{4,2,10},{5,0,2},{5,1,6}});                  c["Xe"]=s;
        s.insert(s.end(), {{4,3,14},{5,2,10},{6,0,2},{6,1,6}});         c["Rn"]=s;
        return c;
    }();
    return cores;
}
int Electrons(const ShellOccupations& s) {int n=0; for (const auto& o : s) n+=o.count; return n;}
} // namespace

int StatedConfiguration::CoreElectrons () const {return Electrons(core);}
int StatedConfiguration::NamedElectrons() const {return Electrons(named);}
ShellOccupations StatedConfiguration::All() const {ShellOccupations a=core; a.insert(a.end(), named.begin(), named.end()); return a;}

Outcome<StatedConfiguration, std::string> ParseConfiguration(const std::string& text)
{
    using R=Outcome<StatedConfiguration, std::string>;
    StatedConfiguration c;
    std::string rest=text;
    const size_t open=rest.find('[');
    if (open!=std::string::npos)
    {
        const size_t close=rest.find(']', open);
        if (close==std::string::npos) return R::Fail("configuration \""+text+"\": an unclosed [ core");
        c.coreSymbol=rest.substr(open+1, close-open-1);
        const auto it=NobleGasCores().find(c.coreSymbol);
        if (it==NobleGasCores().end()) return R::Fail("configuration \""+text+"\": the core ["+c.coreSymbol+"] is not a noble gas (He Ne Ar Kr Xe Rn)");
        c.core=it->second;
        rest=rest.substr(0, open)+" "+rest.substr(close+1);
    }
    std::istringstream is(rest);
    for (std::string tok; is>>tok;)
    {
        size_t i=0;
        while (i<tok.size() && std::isdigit((unsigned char)tok[i])) ++i;
        const char* letter = i<tok.size() ? std::strchr(kLetters, std::tolower((unsigned char)tok[i])) : nullptr;
        if (i==0 || !letter || *letter=='\0') return R::Fail("configuration \""+text+"\": '"+tok+"' is not <n><s|p|d|f>^<count>");
        const int n=std::stoi(tok.substr(0,i)), l=int(letter-kLetters);
        size_t j=i+1;
        if (j<tok.size() && tok[j]=='^') ++j;
        if (j>=tok.size() || !std::all_of(tok.begin()+long(j), tok.end(), [](char ch){return std::isdigit((unsigned char)ch)!=0;}))
            return R::Fail("configuration \""+text+"\": '"+tok+"' has no electron count");
        const int count=std::stoi(tok.substr(j));
        if (n<=l) return R::Fail("configuration \""+text+"\": "+tok+" -- n must exceed l");
        if (count<=0 || count>SubshellCapacity(l))
            return R::Fail("configuration \""+text+"\": "+tok+" -- a "+std::string(1,kLetters[l])+" subshell holds 1 to "+std::to_string(SubshellCapacity(l)));
        for (const auto& o : c.All())
            if (o.n==n && o.l==l) return R::Fail("configuration \""+text+"\": "+std::to_string(n)+kLetters[l]+" is written twice (or is in the core)");
        c.named.push_back({n,l,count});
    }
    if (c.All().empty()) return R::Fail("configuration \""+text+"\": no subshells");
    return R::Ok(c);
}

std::string WhyNotRepresentable(const ShellOccupations& shells, bool fromLowestN)
{
    for (int l=0;l<=LMax;l++)
    {
        std::vector<ShellOccupation> sl;
        for (const auto& o : shells) if (o.l==l) sl.push_back(o);
        if (sl.empty()) continue;
        std::sort(sl.begin(), sl.end(), [](const auto& a, const auto& b){return a.n<b.n;});
        const std::string L(1, kLetters[l]);
        if (fromLowestN && sl.front().n!=l+1)
            return "the lowest occupied "+L+" level is "+std::to_string(sl.front().n)+L+", but an all-electron atom fills "+std::to_string(l+1)+L+" first";
        for (size_t i=0;i+1<sl.size();i++)
        {
            if (sl[i+1].n!=sl[i].n+1)
                return "a hole below an occupied level: "+std::to_string(sl[i].n)+L+" then "+std::to_string(sl[i+1].n)+L
                       +" -- the radial atom fills each l block from its lowest level up";
            if (sl[i].count!=SubshellCapacity(l))
                return std::to_string(sl[i].n)+L+"^"+std::to_string(sl[i].count)+" is open below the occupied "+std::to_string(sl[i+1].n)+L
                       +" -- the radial atom fills each l block from its lowest level up, so only the highest "+L+" subshell may be open";
        }
    }
    return {};
}

std::string ConfigurationText(const ShellOccupations& shells)
{
    ShellOccupations s=shells;
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b){return a.n!=b.n ? a.n<b.n : a.l<b.l;});
    std::string t;
    for (const auto& o : s) t+=(t.empty() ? "" : " ")+std::to_string(o.n)+kLetters[o.l]+"^"+std::to_string(o.count);
    return t;
}

Atom_EC::Atom_EC(const ShellOccupations& shells)
: itsLMax(0), itsLValance(0)
{
    if (const std::string why=WhyNotRepresentable(shells, false); !why.empty())
        throw std::logic_error("Atom_EC: the stated configuration "+ConfigurationText(shells)+" is not representable ("+why+") -- the caller checks "
                               "WhyNotRepresentable first");
    for (const auto& o : shells)
    {
        if (o.l<0 || o.l>LMax) throw std::logic_error("Atom_EC: l="+std::to_string(o.l)+" is past LMax");
        if (o.count==SubshellCapacity(o.l)) itsNs.Nf[o.l]+=o.count;   // closed: core
        else { itsNs.Nv[o.l]=o.count; itsNs.Nu[o.l]=std::min(o.count, SubshellCapacity(o.l)-o.count); }   // the l's one open subshell: Hund
    }
    for (int l=0;l<=LMax;l++)
    {
        itsNs.N[l]=itsNs.Nf[l]+itsNs.Nv[l];
        if (itsNs.N[l]>0) itsLMax=l;
        if (itsNs.Nv[l]>0) itsLValance=l;
    }
    BuildNROccupations();
}

ShellOccupations Atom_EC::Shells(int coreElectrons) const
{
    // The core's shells per l, filled in aufbau order (a pseudo-atom's core is the noble gas below its valence, semicore excluded).
    static const int aufbau[][2]={{1,0},{2,0},{2,1},{3,0},{3,1},{4,0},{3,2},{4,1},{5,0},{4,2},{5,1},{6,0},{4,3},{5,2},{6,1},{7,0}};
    int coreShells[LMax+1]={0,0,0,0};
    int left=coreElectrons;
    for (const auto& sh : aufbau)
    {
        if (left<=0) break;
        left-=SubshellCapacity(sh[1]); coreShells[sh[1]]++;
    }
    if (left!=0) throw std::logic_error("Atom_EC::Shells: a core of "+std::to_string(coreElectrons)+" electrons is not whole aufbau shells");
    ShellOccupations out;
    for (int l=0;l<=LMax;l++)
    {
        int k=itsNs.N[l], n=l+1+coreShells[l];
        while (k>0) { const int c=std::min(k, SubshellCapacity(l)); out.push_back({n,l,c}); k-=c; ++n; }
    }
    return out;
}

// Pseudo-atom: load ONLY the raw valence shells (no core, no full-subshell demotion), so a closed valence
// subshell like Si 3s^2 stays a valence occupation instead of folding into the core.  netCharge != 0 makes
// an ion (see the header): adjust the neutral valence by -netCharge electrons, then recompute the unpaired
// count by Hund's rule on the (adjusted) valence shell.
namespace {
//! The SEMICORE closed shells a pseudopotential of valence \a zion keeps beyond the neutral valence \a nValence, as electron counts per
//! l.  They are the OUTERMOST shells of the noble-gas core (Z - nValence electrons, filled in aufbau order), taken highest n first then
//! highest l -- the energetic order of a core (Ag q19 keeps 4s^2 4p^6, not 3d).  THROWS when \a zion does not peel whole shells.
std::vector<int> SemicoreShells(int Z, int nValence, int zion)
{
    std::vector<int> nf(LMax+1, 0);
    int extra=zion-nValence;
    if (extra==0) return nf;
    if (extra<0) throw std::invalid_argument("PseudoAtom_EC: Z="+std::to_string(Z)+" q="+std::to_string(zion)+" keeps FEWER electrons than "
                                             "the neutral valence ("+std::to_string(nValence)+")");
    static const int aufbau[][2]={{1,0},{2,0},{2,1},{3,0},{3,1},{4,0},{3,2},{4,1},{5,0},{4,2},{5,1},{6,0},{4,3},{5,2},{6,1},{7,0}};
    std::vector<std::pair<int,int>> core;   // (n,l)
    int nCore=0;
    for (const auto& sh : aufbau) { if (nCore>=Z-nValence) break; core.push_back({sh[0],sh[1]}); nCore+=2*(2*sh[1]+1); }
    if (nCore!=Z-nValence) throw std::logic_error("PseudoAtom_EC: the core of Z="+std::to_string(Z)+" is not whole aufbau shells");
    std::sort(core.begin(), core.end(), [](const auto& a, const auto& b){return a.first!=b.first ? a.first>b.first : a.second>b.second;});
    for (const auto& [n,l] : core)
    {
        if (extra==0) break;
        const int g=2*(2*l+1);
        if (g>extra) throw std::invalid_argument("PseudoAtom_EC: Z="+std::to_string(Z)+" q="+std::to_string(zion)+" is not a whole-shell "
                                                 "semicore split (the next core shell n="+std::to_string(n)+" l="+std::to_string(l)+" holds "
                                                 +std::to_string(g)+", "+std::to_string(extra)+" left to place)");
        nf[l]+=g; extra-=g;
    }
    if (extra>0) throw std::invalid_argument("PseudoAtom_EC: q="+std::to_string(zion)+" exceeds Z="+std::to_string(Z));
    return nf;
}
} // namespace

Atom_EC::Atom_EC(int Z, int netCharge, int zion, ValenceOnly_t)
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
    for (int rem = netCharge; rem>0; --rem)          // cation: physical (aufbau-reverse) removal order
    {
        // A transition-metal / rare-earth valence (occupied d/f) holds its s (and p) electrons in the
        // HIGHER principal shell (4s vs 3d), and a cation empties that higher-n shell FIRST: Mn2+ is
        // 3d^5 4s^0 (the S=5/2 MnO/LiMn2O4 ion), NOT 4s^2 3d^3.  A main-group valence (s/p only, same n)
        // removes highest-l first as before (Al2+ -> 3s^1).
        int l = -1;
        if (nv[2]>0 || nv[3]>0)
            { l=0; while (l<=1 && nv[l]<=0) ++l; if (l>1 || nv[l]<=0) l=-1; }   // s then p, if occupied
        if (l<0)
            { l=int(LMax); while (l>=0 && nv[l]<=0) --l; }                      // else highest l down
        assert(l>=0 && "PseudoAtom_EC: cation removes more electrons than the valence holds");
        --nv[l];
    }

    // The SEMICORE the pseudopotential keeps (closed shells below the valence): the ion was made on the valence shells above, so a
    // cation empties 3s before the semicore 2p is touched (Na+ q9 = 2s^2 2p^6).
    int nNeutral=0;
    for (size_t l=0;l<=LMax;l++) nNeutral+=int(vc[l]);
    const std::vector<int> nf=SemicoreShells(Z, nNeutral, zion<0 ? nNeutral : zion);
    for (size_t l=0;l<=LMax;l++)
    {
        itsNs.Nf[l]=nf[l];             //the semicore shells the PP keeps (0 for a valence-only q): the rest of the core is the PP's
        itsNs.Nv[l]=nv[l];
        itsNs.N [l]=nf[l]+nv[l];
        if (itsNs.N[l]>0) itsLMax=l;
        if (nv[l]>0) itsLValance=l;
    }
    // Neutral atom: keep the periodic table's exact unpaired count (handles half-filled-shell exceptions like
    // Cr/Cu).  Ion: recompute by Hund's rule PER SHELL and sum -- n in 2g spin-orbitals gives n unpaired if
    // n<=g else 2g-n -- so a closed-shell ion (F-/Na+/O2-) is 0, Mn2+ d^5 is 5, and a multi-open-shell ion
    // (Mn+ 4s^1 3d^5, the ^7S ground state) counts every open shell (1+5=6).
    int nup;
    if (netCharge==0)
        nup = (int)thePeriodicTable().GetNumUnpairedElectrons(Z);
    else
    {
        nup = 0;
        for (size_t l=0;l<=LMax;l++)
        {
            int g=2*int(l)+1;
            nup += (nv[l]<=g) ? nv[l] : (2*g-nv[l]);
        }
    }
    AssignUnpaired(Z, nup);
}

PseudoAtom_EC::PseudoAtom_EC(int Z, int netCharge, int zion)
: Atom_EC(Z, netCharge, zion, ValenceOnly_t{})
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
   

} // namespace qchem