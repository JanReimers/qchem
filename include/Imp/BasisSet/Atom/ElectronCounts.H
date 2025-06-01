// File:  ElectronCounts.H  Simple structure that store electron configuration counts for each l state.
#ifndef _ElectronCounts_H_
#define _ElectronCounts_H_

const int LMax=3;

struct ElCounts
{
    ElCounts() : N{0,0,0,0}, Nf{0,0,0,0},Nv{0,0,0,0},Nu{0,0,0,0} {};
    int N [LMax+1]; //Total N[l]
    int Nf[LMax+1]; //Full shell N[l]
    int Nv[LMax+1]; //Valance shell N[l]
    int Nu[LMax+1]; //# of unpaired electrons for a given l.
};

enum class Spin;

struct ElCounts_l
{
    int n;  //# valance electrons for a given l
    int nu; //# of un paired electrons (for a given l)
    int GetN(Spin s) const;
};



#endif //_ElectronCounts_H_