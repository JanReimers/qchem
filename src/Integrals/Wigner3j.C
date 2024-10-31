#include "Imp/Integrals/Wigner3j.H"
#include "wignerSymbols/wignerSymbols-cpp.h"
#include <cassert>
#include <iostream>

Wigner3j Wigner3j::theW3j;

using std::cout;
using std::endl;
Wigner3j::Wigner3j()
{
    std::cout << "Initializing Wigner 3j tables LMax=" << LMax << std::endl;
    // Load up data array with markers.
    for (int la=0; la<=LMax; la++)
        for (int l=0; l<=LMax; l++)
            for (int lb=0; lb<=LMax; lb++)
                Data[la][l][lb]=-1.0; //Use this as a marker for un-assigned.
    //
    //  Now make it toally symmetric
    //
    for (int la=0; la<=LMax; la++)
        for (int l=0; l<=LMax; l++)
            for (int lb=0; lb<=LMax; lb++)
                Data[la][l][lb]=0.5*WignerSymbols::wigner3j(la,l,lb,0,0,0)*WignerSymbols::wigner3j(la,l,lb,0,0,0); //Use this as a marker for un-assigned.
          
  
}

using std::cout;
using std::endl;

double Wigner3j::operator()(int la, int lb, int k) const 
{
    assert(la>=0);
    assert(la<=LMax);
    assert(k >=0);
    assert(k <=LMax);
    assert(lb>=0);
    assert(lb<=LMax);
    double ret=Data[la][k][lb];
    if (ret==-1.0)
    {
        cout << "Wigner3j no data for la,k,lb = " << LMax << " " << la  << " " << k << " " << lb << endl;
    }
    assert(ret!=-1.0);
    assert(ret>0.0);
    return 2*ret;
}


