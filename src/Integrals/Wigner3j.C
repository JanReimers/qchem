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

    for (int la=0; la<=LMax; la++)
        for (int lb=0; lb<=LMax; lb++)
            for (int k=0; k<=2*LMax; k++)
                Data[la][lb][k]=WignerSymbols::wigner3j(la,lb,k,0,0,0); //Use this as a marker for un-assigned.
          
  
}

double Wigner3j::operator()(int la, int lb, int k) const 
{
    assert(la>=0);
    assert(la<=LMax);
    assert(k >=0);
    assert(k <=2*LMax);
    assert(lb>=0);
    assert(lb<=LMax);
    return Data[la][lb][k];
}


