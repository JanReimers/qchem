// File: Symmetry/OkmjQN.C  Spherical Spinor Omega_kmj symmetry.

#include "Imp/BasisSet/Atom/kappa/Okmj.H"
#include "Imp/BasisSet/Atom/EC.H"
#include <iostream>
#include <iomanip>
#include <cassert>

using std::cout;
using std::endl;
Omega_k_Sym::Omega_k_Sym(): kappa(0) {};

Omega_k_Sym::Omega_k_Sym(int _kappa) : kappa(_kappa) 
{
    assert(abs(kappa)<10);
};

size_t Omega_k_Sym::SequenceIndex() const //Used for op<
 {
    assert(abs(kappa)<=LMax+1);
    return kappa+LMax+1;
 }

bool Omega_k_Sym::MatchType(const Symmetry& b) const
{
    return dynamic_cast<const Omega_k_Sym*>(&b)!=0;
}

bool Omega_k_Sym::Match(const Symmetry& qn) const
{
    const Omega_k_Sym* oqn = dynamic_cast<const Omega_k_Sym*>(&qn);
    assert(oqn);
    return kappa==oqn->kappa;
}

int Omega_k_Sym::GetDegeneracy() const
{
    return Getj()+0.5; //(2j+1)/2 degeneracy for one spin state.
}

ElCounts_l Omega_k_Sym::GetN(const ElCounts& ec) const
{
    int l=GetL();
    int nl=ec.N[l];
    if (ec.Nv[l]==0) return ElCounts_l{nl,0};//::make_pair(nl,0);
    assert(nl!=0);
    // Handle partial shells
    int nlu=1; //# unpaired in shell l. 
    int NUnpaired=0;
    for (auto nu:ec.Nu) NUnpaired+=nu;
    if (l==1) // p is partial.
    {
        assert(ec.Nv[2]==0); //No partial D orbital
        nlu=NUnpaired-ec.Nv[0];
    }
    else if (l==2) // d is partial.
    {
        assert(ec.Nv[1]==0); //p better be full
        if (ec.Nv[l]>1) nlu=NUnpaired-ec.Nv[0];            
    }
    else if(l==3) // f is partial.
    {
        
        assert(ec.Nv[0]==0); //If f is Partial s must be full.
        assert(ec.Nv[1]==0); //If f is Partial p must be full.
        nlu=NUnpaired-ec.Nv[2];
        assert(nlu>=0);
    }
    return ElCounts_l{nl,nlu};// std::make_pair(nl,nlu);
}

extern std::string SPDFG[];
std::string j2s[]={"1/2","3/2","5/2","7/2","9/2"};


std::ostream& Omega_k_Sym::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        int jindex=Getj()-0.5;
        os << SPDFG[GetL()] << j2s[jindex] << " kappa=" << std::setw(2) << kappa << " ";
        
    }
    return os;
}

Angular_Sym* Omega_k_Sym::Clone() const
{
    return new Omega_k_Sym(*this);
}



Omega_kmj_Sym::Omega_kmj_Sym(): kappa(0), mj(0) {};

Omega_kmj_Sym::Omega_kmj_Sym(int _kappa, double _mj) : kappa(_kappa), mj(_mj) {};

size_t Omega_kmj_Sym::SequenceIndex() const //Used for op<
 {
    assert(abs(kappa)<=LMax+1);
    return (mj+Getj())*(2*LMax+3)+(kappa+LMax+1);
 }

bool Omega_kmj_Sym::MatchType(const Symmetry& b) const
{
    return dynamic_cast<const Omega_kmj_Sym*>(&b)!=0;
}

bool Omega_kmj_Sym::Match(const Symmetry& qn) const
{
    const Omega_kmj_Sym* oqn = dynamic_cast<const Omega_kmj_Sym*>(&qn);
    assert(oqn);
    return kappa==oqn->kappa && mj==oqn->mj;
}

int Omega_kmj_Sym::GetDegeneracy() const
{
    return 1;
}

//
// Get the electron config assuming all mj states for a given k are degenerate.
//
std::pair<int,int> Omega_kmj_Sym::GetNk(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const
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



std::pair<int,int> Omega_kmj_Sym::GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const
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

ElCounts_l Omega_kmj_Sym::GetN(const ElCounts& ec) const
{
    //assert(itsL<=LMax);
    int NUnpaired=0;
    for (auto nu:ec.Nu) NUnpaired+=nu;
    int nl,nlu;
    std::tie(nl,nlu)=GetNk(ec.N,ec.Nv,NUnpaired);
    assert((nl+nlu)%2==0);
    double j=Getj();
    int g=2*j+1;
    int l=GetL();
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
        int ml=Getml();
        nlv=nlmv[ml+l];
        
        nlu=nlmu[ml+l];
        //cout << "(" << " " << nlmv << " " << nlmu << ") ";
    }
    assert(nlv%2==nlu);
    return ElCounts_l{nlc+nlv,nlu};//::make_pair(nlc+nlv,nlu);
}

std::ostream& Omega_kmj_Sym::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        int jindex=Getj()-0.5;
        os << SPDFG[GetL()] << j2s[jindex] << " kappa=" << std::setw(2) << kappa << " mj=" << std::setw(4) << std::setprecision(1) << mj << " ";
        
    }
    return os;
}

std::istream& Omega_kmj_Sym::Read (std::istream& is)
{
    return is;
}

Angular_Sym* Omega_kmj_Sym::Clone() const
{
    return new Omega_kmj_Sym(*this);
}

