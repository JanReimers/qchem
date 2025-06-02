// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/BasisSet/Atom/ml/Ylm.H"
#include "Imp/BasisSet/Atom/EC.H"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <Spin.H>

using std::cout;
using std::endl;

const int LMAX=4;
Ylm_Sym::Ylm_Sym(): Yl_Sym(0) {};


Ylm_Sym::Ylm_Sym(int l, const std::vector<int>& _ml) 
: Yl_Sym(l),  ml(_ml) 
{
    assert(ml.size()>0);
};

size_t Ylm_Sym::SequenceIndex() const //Used for op<
 {
    int mmax=*std::max_element(ml.begin(),ml.end());
    return mmax+itsL+itsL*(2*LMAX+1);
 }

bool Ylm_Sym::MatchType(const Symmetry& b) const
{
    return dynamic_cast<const Ylm_Sym*>(&b)!=0;
}

bool Ylm_Sym::Match(const Symmetry& qn) const
{
    const Ylm_Sym* yqn = dynamic_cast<const Ylm_Sym*>(&qn);
    assert(yqn);
    bool meq=ml.size()==yqn->ml.size() && std::equal(ml.begin(), ml.end(), yqn->ml.begin());
    return itsL==yqn->itsL && meq;
}

int Ylm_Sym::GetDegeneracy() const
{
    return ml.size(); 
}

//  This function is non trivial  we need to count how many total, and how many unpaired electrons
//  there are for a given l and list of ml values.
//  We do this by dumping elkectrons into valance and un paired buckets.
ElCounts_l Ylm_Sym::GetN(const ElCounts& ec) const
{
    ec.DebugCheck();
    assert(itsL<=LMax);
    int l=itsL;

    int nup=ec.GetNv(l,Spin::Up);   //How many spin up valance electrons.
    int ndn=ec.GetNv(l,Spin::Down); //How many spin dn valance electrons.
    //  Fill all up and down occupied m states.
    int nlmups[2*LMax+1]={0,0,0,0,0,0,0};
    int nlmdns[2*LMax+1]={0,0,0,0,0,0,0};
    for (int m=-l;m<=l && nup>0;m++,nup--) nlmups[m+l]++;
    for (int m=-l;m<=l && ndn>0;m++,ndn--) nlmdns[m+l]++;
    // Now count how many of these occupied states are in the ml list.
    int nlmup=0,nlmdn=0;
    for (auto m:ml)
    {
        nlmup+=nlmups[m+l];
        nlmdn+=nlmdns[m+l];
    }
    assert(nlmup>=nlmdn);
    int nlmv=nlmup+nlmdn; //Total valance
    int nlmu=nlmup-nlmdn; //Total unpaired.

    int g=2*l+1; //Degeneracy.
    int nlmc=ec.Nf[l]/g*ml.size();  //How many m states in the core.

    return ElCounts_l{nlmc+nlmv,nlmu};//// Should be total core+valance 
}


extern std::string SPDFG[];

inline int width(int m) {return m<0 ? 2 : 1;}
std::ostream& Ylm_Sym::Write(std::ostream& os) const
{
    os << SPDFG[itsL] << " ";
    int n_extra=22;
    if (ml.size()<2*(size_t)itsL+1)
    {
        n_extra-=4;
        os << "[ ";
        for (auto im:ml)
        {
            int w=width(im);
            os << std::setw(w) << im << " ";
            n_extra-=(w+1);
        }
        os << "] ";

    }
    for (int i=0;i<n_extra;i++) os << " ";
    os << "|";
    return os;
}

Angular_Sym* Ylm_Sym::Clone() const
{
    return new Ylm_Sym(*this);
}

