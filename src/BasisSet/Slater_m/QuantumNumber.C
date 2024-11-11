// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/BasisSet/Slater_m/QuantumNumber.H"
#include <iostream>
#include <cassert>

using std::cout;
using std::endl;

YlmQN::YlmQN(): SphericalSymmetryQN(0),m(0) {};

YlmQN::YlmQN(int _l, int _m) : SphericalSymmetryQN(_l), m(_m) {};

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

std::pair<int,int> YlmQN::GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const
{
    int nl,nlu;
    std::tie(nl,nlu)=SphericalSymmetryQN::GetN(N,Nv,NUnpaired);
    assert((nl+nlu)%2==0);
    int l=itsL;
    int g=2*l+1;
    int nlc=N[l]-Nv[l];
    assert(nlc%(2*g)==0);
    nlc/=g;
    int nlv=Nv[l];

    if (l==1)
    {
        int nlmv[3]={0,0,0};
        int nlmu[3]={0,0,0};
        for (int m1=-1;m1<=l&&nlv>0;m1++)
        {
            nlmv[m1+l]++;
            nlmu[m1+l]++;
            nlv--;
        }
        for (int m1=-1;m1<=l&&nlv>0;m1++)
        {
            nlmv[m1+l]++;
            nlmu[m1+l]--;
            nlv--;
        }
        nlv=nlmv[m+l];
        nlu=nlmu[m+l];
        //cout << "(" << " " << nlmv << " " << nlmu << ") ";
    }
    if (l==2)
    {
        int nlmv[5]={0,0,0,0,0};
        int nlmu[5]={0,0,0,0,0};
        for (int m1=-l;m1<=l&&nlv>0;m1++)
        {
            nlmv[m1+l]++;
            nlmu[m1+l]++;
            nlv--;
        }
        for (int m1=-l;m1<=l&&nlv>0;m1++)
        {
            nlmv[m1+l]++;
            nlmu[m1+l]--;
            nlv--;
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
        os << SPDFG[itsL] << "_" << m << " ";
    return os;
}

std::istream& YlmQN::Read (std::istream& is)
{
    return is;
}

QuantumNumber* YlmQN::Clone() const
{
    return new YlmQN(*this);
}

