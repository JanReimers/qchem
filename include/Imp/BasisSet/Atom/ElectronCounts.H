// File:  ElectronCounts.H  Simple structure that store electron configuration counts for each l state.
#ifndef _ElectronCounts_H_
#define _ElectronCounts_H_

const int LMax=3;

struct ElCounts
{
    ElCounts() : N{0,0,0,0}, Nf{0,0,0,0},Nv{0,0,0,0},Nu{0,0,0,0} {};
    void DebugCheck() const; //Check self consistency

    int N [LMax+1]; //Total N[l]
    int Nf[LMax+1]; //Full shell or core, N[l]
    int Nv[LMax+1]; //Valance shell N[l].  He is considered to have no valance electrons ... all are core.
    int Nu[LMax+1]; //# of unpaired electrons for a given l.
};

enum class Spin;

struct ElCounts_l
{
    int N;  //# Total electrons for a given l
    int Nu; //# of un paired electrons (for a given l)
    int GetN(Spin s) const;
};



#endif //_ElectronCounts_H_