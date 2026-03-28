// File: AngularIntegrals.C  Handle the angular part of 2-electron ERIs
module;
#include <blaze/Math.h>
export module qchem.BasisSet.Atom.Internal.AngularIntegrals;

export namespace AngularIntegrals
{
    typedef blaze::StaticVector<double,11> rvec11_t; //la+lc+1<=11 support up to i orbitals.  Good luck finding a stable nucleus!

    double Exchange(int k,int la,int lb);
    double Coulomb (int k,int la,int lc,int ma,int mc);            
    double Exchange(int k,int la,int lb,int ma,int mb);            

    rvec11_t Coulomb (int la,int lc);            
    rvec11_t Exchange(int la,int lb);            
    rvec11_t Coulomb (int la,int lc,int ma,int mc);            
    rvec11_t Exchange(int la,int lb,int ma,int mb);            

};
