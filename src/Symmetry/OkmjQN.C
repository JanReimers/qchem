// File: Symmetry/OkmjQN.C  Spherical Spinor Omega_kmj symmetry.

#include "Imp/Symmetry/OkmjQN.H"
#include "Imp/Symmetry/AtomQN.H"
#include <iostream>
#include <iomanip>
#include <cassert>

using std::cout;
using std::endl;
Omega_kQN::Omega_kQN(): kappa(0) {};

Omega_kQN::Omega_kQN(int _kappa) : kappa(_kappa) {};

bool Omega_kQN::Match(const QuantumNumber& qn) const
{
    const Omega_kQN* oqn = dynamic_cast<const Omega_kQN*>(&qn);
    assert(oqn);
    return kappa==oqn->kappa;
}

int Omega_kQN::GetDegeneracy() const
{
    return Getj()+0.5; //(2j+1)/2 degeneracy for one spin state.
}

QuantumNumber* Omega_kQN::AddPrincipleQN(int index) const
{
    return new AtomQN(index,*this);
}

std::pair<int,int> Omega_kQN::GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const
{
    int l=GetL();
    int nl=N[l];
    if (Nv[l]==0) return std::make_pair(nl,0);
    assert(nl!=0);
    // Handle partial shells
    int nlu=1; //# unpaired in shell l. 
    if (l==1) // p is partial.
    {
        assert(Nv[2]==0); //No partial D orbital
        nlu=NUnpaired-Nv[0];
    }
    else if (l==2) // d is partial.
    {
        assert(Nv[1]==0); //p better be full
        if (Nv[l]>1) nlu=NUnpaired-Nv[0];            
    }
    else if(l==3) // f is partial.
    {
        
        assert(Nv[0]==0); //If f is Partial s must be full.
        assert(Nv[1]==0); //If f is Partial p must be full.
        nlu=NUnpaired-Nv[2];
        assert(nlu>=0);
    }
    return std::make_pair(nl,nlu);
}

extern std::string SPDFG[];
std::string j2s[]={"1/2","3/2","5/2","7/2","9/2"};


std::ostream& Omega_kQN::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        int jindex=Getj()-0.5;
        os << SPDFG[GetL()] << j2s[jindex] << " kappa=" << std::setw(2) << kappa << " ";
        
    }
    return os;
}

std::istream& Omega_kQN::Read (std::istream& is)
{
    return is;
}

AngularQN* Omega_kQN::Clone() const
{
    return new Omega_kQN(*this);
}







Omega_kmjQN::Omega_kmjQN(): kappa(0), mj(0) {};

Omega_kmjQN::Omega_kmjQN(int _kappa, double _mj) : kappa(_kappa), mj(_mj) {};

bool Omega_kmjQN::Match(const QuantumNumber& qn) const
{
    const Omega_kmjQN* oqn = dynamic_cast<const Omega_kmjQN*>(&qn);
    assert(oqn);
    return kappa==oqn->kappa && mj==oqn->mj;
}

int Omega_kmjQN::GetDegeneracy() const
{
    return 1;
}

QuantumNumber* Omega_kmjQN::AddPrincipleQN(int index) const
{
    return new AtomQN(index,*this);
}

//
// Get the electron config assuming all mj states for a given k are degenerate.
//
std::pair<int,int> Omega_kmjQN::GetNk(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const
{
    int l=GetL();
    int nl=N[l];
    if (Nv[l]==0) return std::make_pair(nl,0);
    assert(nl!=0);
    // Handle partial shells
    int nlu=1; //# unpaired in shell l. 
    if (l==1) // p is partial.
    {
        assert(Nv[2]==0); //No partial D orbital
        nlu=NUnpaired-Nv[0];
    }
    else if (l==2) // d is partial.
    {
        assert(Nv[1]==0); //p better be full
        if (Nv[l]>1) nlu=NUnpaired-Nv[0];            
    }
    else if(l==3) // f is partial.
    {
        
        assert(Nv[0]==0); //If f is Partial s must be full.
        assert(Nv[1]==0); //If f is Partial p must be full.
        nlu=NUnpaired-Nv[2];
        assert(nlu>=0);
    }
    return std::make_pair(nl,nlu);
}


std::pair<int,int> Omega_kmjQN::GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const
{
    //assert(itsL<=LMax);
    int nl,nlu;
    std::tie(nl,nlu)=GetNk(N,Nv,NUnpaired);
    assert((nl+nlu)%2==0);
    double j=Getj();
    int g=2*j+1;
    int l=GetL();
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
        int ml=Getml();
        nlv=nlmv[ml+l];
        
        nlu=nlmu[ml+l];
        //cout << "(" << " " << nlmv << " " << nlmu << ") ";
    }
    assert(nlv%2==nlu);
    return std::make_pair(nlc+nlv,nlu);
}



std::ostream& Omega_kmjQN::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        int jindex=Getj()-0.5;
        os << SPDFG[GetL()] << j2s[jindex] << " kappa=" << std::setw(2) << kappa << " mj=" << std::setw(4) << std::setprecision(1) << mj << " ";
        
    }
    return os;
}

std::istream& Omega_kmjQN::Read (std::istream& is)
{
    return is;
}

AngularQN* Omega_kmjQN::Clone() const
{
    return new Omega_kmjQN(*this);
}

