// File: testherm.cpp

#include "Hermite/Hermite2.H"
#include "dft/Polarization.H"
#include <iostream.h>

int main()
{
  StreamableObject::SetOutputMode(StreamableObject::pretty);
  {
    Hermite2 d(1.0,RVec(-1,0,0),RVec(1,0,0),1,2);
    cout << d(Polarization(0,0,0),Polarization(0,0,1),Polarization(2,0,0)) << endl << d;
  }
  {
    Hermite2 d(1.0,RVec(-1,0,0),RVec(1,0,0),2,3);
    cout << d(Polarization(0,0,0),Polarization(0,0,1),Polarization(2,0,0)) << endl << d;
  }
 
  
  return 0;
}
