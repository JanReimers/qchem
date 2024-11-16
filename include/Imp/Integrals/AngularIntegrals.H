// File: AngularIntegrals.H  Handle the angular part of 2-electron ERIs
#ifndef _AngularIntegrals_H_
#define _AngularIntegrals_H_

template <class T> class Vector;

namespace AngularIntegrals
{
    double Exchange(int k,int la,int lb);
    double Coulomb (int k,int la,int lc,int ma,int mc);            
    double Exchange(int k,int la,int lb,int ma,int mb);            

    typedef Vector<double> RVec;
    RVec Coulomb (int la,int lc);            
    RVec Exchange(int la,int lb);            
    RVec Coulomb (int la,int lc,int ma,int mc);            
    RVec Exchange(int la,int lb,int ma,int mb);            

};


#endif // _AngularIntegrals_H_
