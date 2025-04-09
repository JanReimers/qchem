// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/Symmetry/Ylm.H"
#include "Imp/WaveFunction/Atom_EC.H"
#include <iostream>
#include <iomanip>
#include <cassert>

using std::cout;
using std::endl;

const int LMAX=4;
Ylm_Sym::Ylm_Sym(): Yl_Sym(0),m(0) {};

Ylm_Sym::Ylm_Sym(int _l, int _m) : Yl_Sym(_l), m(_m) {};

size_t Ylm_Sym::SequenceIndex() const //Used for op<
 {
    return m+itsL+itsL*(2*LMAX+1);
 }

bool Ylm_Sym::MatchType(const Symmetry& b) const
{
    return dynamic_cast<const Ylm_Sym*>(&b)!=0;
}

bool Ylm_Sym::Match(const Symmetry& qn) const
{
    const Ylm_Sym* yqn = dynamic_cast<const Ylm_Sym*>(&qn);
    assert(yqn);
    return itsL==yqn->itsL && m==yqn->m;
}

int Ylm_Sym::GetDegeneracy() const
{
    return 1; 
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
    nlc/=g;
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
        nlv=nlmv[m+l];
        
        nlu=nlmu[m+l];
        //cout << "(" << " " << nlmv << " " << nlmu << ") ";
    }
    assert(nlv%2==nlu);
    return ElCounts_l{nlc+nlv,nlu};//::make_pair(nlc+nlv,nlu);
}


extern std::string SPDFG[];

std::ostream& Ylm_Sym::Write(std::ostream& os) const
{
    return os << SPDFG[itsL] << " " << std::setw(2) << m << " ";
}

Angular_Sym* Ylm_Sym::Clone() const
{
    return new Ylm_Sym(*this);
}

