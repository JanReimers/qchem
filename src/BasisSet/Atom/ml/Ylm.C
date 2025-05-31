// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/BasisSet/Atom/ml/Ylm.H"
#include "Imp/BasisSet/Atom/EC.H"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>

using std::cout;
using std::endl;

const int LMAX=4;
Ylm_Sym::Ylm_Sym(): Yl_Sym(0) {};


Ylm_Sym::Ylm_Sym(int l, const std::vector<int>& _ml) 
: Yl_Sym(l),  ml(_ml) 
{


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

ElCounts_l Ylm_Sym::GetN(const ElCounts& ec) const
{
    assert(itsL<=LMax);
    ElCounts_l ecl=Yl_Sym::GetN(ec);
    int l=itsL;
    int g=2*l+1;

    int nlv=ec.Nv[l]; //#valance for this l and all m
    int nlu=ecl.nu; //# unpaired for this l and all m
    int nlc=ec.N[l]-ec.Nv[l]; //#core for this l and all m
    assert((ecl.n+ecl.nu)%2==0);
    assert(nlc%(2*g)==0);

 

    
    int nlmvs[2*LMax+1]={0,0,0,0,0,0,0};
    int nlmus[2*LMax+1]={0,0,0,0,0,0,0};
    int nlmv=0,nlmu=0;
   
        //cout << "Start v,u=" << nlv << " " << nlu << endl;

    bool less_than_half = nlv<=g;
    for (int m1=-l;m1<=l && nlv>0 && nlu>=0;m1++)
    {
        nlmvs[m1+l]++;
        nlmus[m1+l]++;
        nlv--;
        if (less_than_half) nlu--;
        //cout << "Up v,u=" << nlv << " " << nlu << endl;
    }
    for (int m1=-l;m1<=l&&nlv>0;m1++)
    {
        nlmvs[m1+l]++;
        nlmus[m1+l]--;
        nlv--;
        if (!less_than_half) nlu--;
        //cout << "Down v,u=" << nlv << " " << nlu << endl;
    }

    for (auto im:ml)
    {
        nlmv+=nlmvs[im+l];
        nlmu+=nlmus[im+l];
    }
        //cout << "(" << " " << nlmv << " " << nlmu << ") ";
    int nlmc=nlc/g*ml.size();

    return ElCounts_l{nlmc+nlmv,nlmu};//::make_pair(nlc+nlv,nlu);
}


extern std::string SPDFG[];

std::ostream& Ylm_Sym::Write(std::ostream& os) const
{
    os << SPDFG[itsL] << " ";
    for (auto im:ml)
        os << std::setw(2) << im << " ";

    return os;
}

Angular_Sym* Ylm_Sym::Clone() const
{
    return new Ylm_Sym(*this);
}

