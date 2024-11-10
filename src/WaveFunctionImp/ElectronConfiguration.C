// File: ElectronConfiguration.C

#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/Misc/PeriodicTable.H"
#include <Spin.H>
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include <cassert>
#include <iostream>
#include <initializer_list>

using std::cout;
using std::endl;

PeriodicTable pt;

const int AtomElectronConfiguration::FullShells[Nshell][LMax+2]=
{
    {0  ,0,0,0,0},
    {2  ,1,0,0,0},
    {10 ,2,1,0,0},
    {18 ,3,2,0,0},
    {36 ,4,3,1,0},
    {54 ,5,4,2,0},
    {86 ,6,5,3,1},
    {118,7,6,4,2}
};

AtomElectronConfiguration::AtomElectronConfiguration(int Z)
: N{0,0,0,0}
, NUnpaired(pt.GetNumUnpairedElectrons(Z))
{
    assert(Z>0);
    assert(Z<=N_Elements);
    int* valance_configuration=pt.GetValanceConfiguration(Z);
    int ns=0;
    for (;ns<Nshell;ns++)
        if (FullShells[ns][0]>Z) break;
    ns--;
//    for (auto l:{0,1,2,3}) cout << valance_configuration[l] << ",";
//    cout << endl;
//    for (auto l:{0,1,2,3}) cout << FullShells[ns][l+1] << ",";
//    cout << endl;

    for (int l=0;l<=LMax;l++) N[l]=FullShells[ns][l+1]*2*(2*l+1)+valance_configuration[l];
//    for (auto l:{0,1,2,3}) cout << N[l] << ",";
//    cout << endl;
}

int AtomElectronConfiguration::GetN() const
{
    int ne=0;
    for (auto n:N) ne+=n;
    return ne;
}

int AtomElectronConfiguration::GetN(const Spin& s) const
{
    int ne=GetN();
    assert((ne+NUnpaired)%2==0);
    assert(s!=Spin::None);
    return s==Spin::Up ? (ne+NUnpaired)/2 : (ne-NUnpaired)/2;
}

int AtomElectronConfiguration::GetN(const QuantumNumber& qn) const
{
    const SphericalSymmetryQN& sqn=dynamic_cast<const SphericalSymmetryQN&>(qn);
    return N[sqn.GetL()];    
}
int AtomElectronConfiguration::GetN(const QuantumNumber& qn, const Spin& s) const
{
    const SphericalSymmetryQN& sqn=dynamic_cast<const SphericalSymmetryQN&>(qn);
    int l=sqn.GetL();
    int nl=N[l];
    if (nl==0) return nl;
    if (l==0)
    { //s
        if (nl%2==0) 
            return nl/2;
        else
            return s==Spin::Up ? (nl+1)/2 : (nl-1)/2;
    }
    else if (l==1)
    {
        if (nl%6==0) return nl/2;
        assert(N[2]%10==0); //No partial D orbital
        int NpUnpaired= N[0]%2==0 ? NUnpaired : NUnpaired-1;
        assert((nl+NpUnpaired)%2==0);
        return s==Spin::Up ? (nl+NpUnpaired)/2 : (nl-NpUnpaired)/2;            
    }
    else if (l==2)
    {
        if (nl%10==0) return nl/2;
        assert(N[1]%6==0); //p better be full
        int ndv=nl%10;
        int nsv=N[0]%2;
//        cout << "nsv,ndv,nfv = " << nsv << " " << ndv << " " << nfv << endl;
        if (nsv==1)
        {
            int NdUnpaired= NUnpaired-1; //One of the unpaired is s?
            return s==Spin::Up ? (nl+NdUnpaired)/2 : (nl-NdUnpaired)/2;            
        }
        if (ndv==1)
        {
            return s==Spin::Up ? (nl+1)/2 : (nl-1)/2;
        }
        if (ndv>1)
            return s==Spin::Up ? (nl+NUnpaired)/2 : (nl-NUnpaired)/2;   
    }
    else if(l==3)
    {
        if (nl%14==0) return nl/2;
        // f is partial.
        assert(N[0]%2==0); //If f is Partial s must be full.
        assert(N[1]%6==0); //If f is Partial p must be full.
        if (N[2]%10==0) //d is full so all unpaired must be f.
        {
            assert((nl+NUnpaired)%2==0);
            return s==Spin::Up ? (nl+NUnpaired)/2 : (nl-NUnpaired)/2;            
        }
        else
        {
            assert(N[2]%10==1); //Only case is one unpaired d electron.
            int NfUnpaired=NUnpaired-1;
            assert((nl+NfUnpaired)%2==0);
            return s==Spin::Up ? (nl+NfUnpaired)/2 : (nl-NfUnpaired)/2;            
        }
    }
    assert(false);
    return 0;
}
    
