// File: BasisSet/Atom/Evaluators/Internal/AngularIntegrals.C  Handle the angular part of 2-electron ERIs
module;
export module qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;
export import qchem.Types;

export namespace AngularIntegrals
{
    double Exchange(int k,int la,int lb);
    double Direct  (int k,int la,int lc,int ma,int mc);            
    double Exchange(int k,int la,int lb,int ma,int mb);            

    rvec11_t Direct  (int la,int lc);            
    rvec11_t Exchange(int la,int lb);            
    rvec11_t Direct  (int la,int lc,int ma,int mc);            
    rvec11_t Exchange(int la,int lb,int ma,int mb);            

};
