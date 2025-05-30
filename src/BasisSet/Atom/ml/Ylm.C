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
Ylm_Sym::Ylm_Sym(): Yl_Sym(0),m(0) {};

Ylm_Sym::Ylm_Sym(int l, int _m) : Yl_Sym(l), m(_m) {};

Ylm_Sym::Ylm_Sym(int l, const std::vector<int>& _ml) : Yl_Sym(l), m(_ml.front()), ml(_ml) {};

size_t Ylm_Sym::SequenceIndex() const //Used for op<
 {
    int mmax=m;
    if (ml.size()>0)
        mmax=*std::max_element(ml.begin(),ml.end());
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
    bool meq =  m==yqn->m;
    if (ml.size()>0)
        meq=ml.size()==yqn->ml.size() && std::equal(ml.begin(), ml.end(), yqn->ml.begin());
    return itsL==yqn->itsL && meq;
}

int Ylm_Sym::GetDegeneracy() const
{
    return ml.size()==0 ? 1 : ml.size(); 
}

ElCounts_l Ylm_Sym::GetN(const ElCounts& ec) const
{
    assert(itsL<=LMax);
    ElCounts_l ecl=Yl_Sym::GetN(ec);
    int nlu=ecl.nu;
    assert((ecl.n+ecl.nu)%2==0);
    int l=itsL;
    int g=2*l+1;
    int nlc=ec.N[l]-ec.Nv[l];
    assert(nlc%(2*g)==0);
    // if (ml.size()==0)
        nlc/=g;
    if (ml.size()>0)
        nlc*=ml.size();

    int nlv=ec.Nv[l];
    
    int nlmv[2*LMax+1]={0,0,0,0,0,0,0};
    int nlmu[2*LMax+1]={0,0,0,0,0,0,0};

    if (l>0)
    {
        //cout << "Start v,u=" << nlv << " " << nlu << endl;

        bool less_than_half = nlv<=g;
        for (int m1=-l;m1<=l&&nlv>0&&nlu>=0;m1++)
        {
            nlmv[m1+l]++;
            nlmu[m1+l]++;
            nlv--;
            if (less_than_half) nlu--;
            //cout << "Up v,u=" << nlv << " " << nlu << endl;
        }
        for (int m1=-l;m1<=l&&nlv>0;m1++)
        {
            nlmv[m1+l]++;
            nlmu[m1+l]--;
            nlv--;
            if (!less_than_half) nlu--;
            //cout << "Down v,u=" << nlv << " " << nlu << endl;
        }
        if (ml.size()==0)
        {
            nlv=nlmv[m+l];
            nlu=nlmu[m+l];
            assert(nlv%2==nlu);
        }
        else
        {
            nlv=0;
            nlu=0;
            for (auto im:ml)
            {
                nlv+=nlmv[im+l];
                nlu+=nlmu[im+l];
            }
        }
        //cout << "(" << " " << nlmv << " " << nlmu << ") ";
    }
    return ElCounts_l{nlc+nlv,nlu};//::make_pair(nlc+nlv,nlu);
}


extern std::string SPDFG[];

std::ostream& Ylm_Sym::Write(std::ostream& os) const
{
    os << SPDFG[itsL] << " ";
    if (ml.size()==0)
        os << std::setw(2) << m << " ";
    else
        for (auto im:ml)
            os << std::setw(2) << im << " ";

    return os;
}

Angular_Sym* Ylm_Sym::Clone() const
{
    return new Ylm_Sym(*this);
}

