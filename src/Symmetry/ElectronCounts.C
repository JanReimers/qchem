// File:  ElectronCounts.C  Simple structure that store electron configuration counts for each l state.
export module qchem.Symmetry.ElectronCounts;
export import qchem.Symmetry.Spin;

export const int LMax=3;

export struct ElCounts
{
    ElCounts() : N{0,0,0,0}, Nf{0,0,0,0},Nv{0,0,0,0},Nu{0,0,0,0} {};
    void DebugCheck() const; //Check self consistency
    int GetNv(int l, Spin s) const;

    int N [LMax+1]; //Total N[l]
    int Nf[LMax+1]; //Full shell or core, N[l]
    int Nv[LMax+1]; //Valance shell N[l].  He is considered to have no valance electrons ... all are core.
    int Nu[LMax+1]; //# of unpaired electrons for a given l.
};


export struct ElCounts_l
{
    int N;  //# Total electrons for a given l
    int Nu; //# of un paired electrons (for a given l)
    int GetN(Spin s) const;
};
