// File: AngularIntegrals.C  Handle the angular part of 2-electron ERIs
export module qchem.BasisSet.Atom.Internal.AngularIntegrals;
export import qchem.Types;

export namespace AngularIntegrals
{
    double Exchange(int k,int la,int lb);
    double Coulomb (int k,int la,int lc,int ma,int mc);            
    double Exchange(int k,int la,int lb,int ma,int mb);            

    RVec Coulomb (int la,int lc);            
    RVec Exchange(int la,int lb);            
    RVec Coulomb (int la,int lc,int ma,int mc);            
    RVec Exchange(int la,int lb,int ma,int mb);            

};
