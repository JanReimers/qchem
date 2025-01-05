// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/Symmetry/YlmQN.H"
#include "Imp/Symmetry/AtomQN.H"
#include <iostream>
#include <iomanip>
#include <cassert>

using std::cout;
using std::endl;

YlmQN::YlmQN(): YlQN(0),m(0) {};

YlmQN::YlmQN(int _l, int _m) : YlQN(_l), m(_m) {};

bool YlmQN::Match(const QuantumNumber& qn) const
{
    const YlmQN* yqn = dynamic_cast<const YlmQN*>(&qn);
    assert(yqn);
    return itsL==yqn->itsL && m==yqn->m;
}

int YlmQN::GetDegeneracy() const
{
    return 1;
}

QuantumNumber* YlmQN::AddPrincipleQN(int index) const
{
    return new AtomQN(index,*this);
}

std::pair<int,int> YlmQN::GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const
{
    assert(itsL<=LMax);
    int nl,nlu;
    std::tie(nl,nlu)=YlQN::GetN(N,Nv,NUnpaired);
    assert((nl+nlu)%2==0);
    int l=itsL;
    int g=2*l+1;
    int nlc=N[l]-Nv[l];
    assert(nlc%(2*g)==0);
    nlc/=g;
    int nlv=Nv[l];
    
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
    return std::make_pair(nlc+nlv,nlu);
}


extern std::string SPDFG[];

std::ostream& YlmQN::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << SPDFG[itsL] << " " << std::setw(2) << m << " ";
    return os;
}

std::istream& YlmQN::Read (std::istream& is)
{
    return is;
}

AngularQN* YlmQN::Clone() const
{
    return new YlmQN(*this);
}

