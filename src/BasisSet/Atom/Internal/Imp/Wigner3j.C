module;
#include "wignerSymbols/wignerSymbols-cpp.h"
#include <cassert>
#include <iostream>
module qchem.BasisSet.Atom.Internal.Wigner3j;

Wigner3j Wigner3j::w3j;

using std::cout;
using std::endl;

Wigner3j::Wigner3j()
{
    std::cout << "Initializing Wigner 3j tables LMax=" << LMax << std::endl;

    for (int la=0; la<=LMax; la++)
        for (int lb=0; lb<=LMax; lb++)
            for (int k=0; k<=2*LMax; k++)
            {
                Data[la][lb][k]=WignerSymbols::wigner3j(la,lb,k,0,0,0);
                for (int ma=-la;ma<=la;ma++)
                    for (int mb=-lb;mb<=lb;mb++)
                        Data_m[la][lb][k][ma+LMax][mb+LMax]=WignerSymbols::wigner3j(la,lb,k,ma,mb,-ma-mb);
            }
          
  
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

double Wigner3j::operator()(int la, int lb, int k,int ma, int mb) const 
{
    assert(la>=0);
    assert(la<=LMax);
    assert(k >=0);
    assert(k <=2*LMax);
    assert(lb>=0);
    assert(lb<=LMax);
    assert(ma>=-la);
    assert(ma<= la);
    assert(mb>=-lb);
    assert(mb<= lb);
    
    return Data_m[la][lb][k][ma+LMax][mb+LMax];
}

