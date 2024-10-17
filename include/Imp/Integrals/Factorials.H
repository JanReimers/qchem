// File: Factorials.H
#ifndef _Factorial_H_
#define _Factorial_H_

namespace qchem {
 const int NMax=19; //Should enable up to f orbitals.
 extern double DFact[]; //Double factorials 1,3,3*5,3*5*7 etc. lookup table.
 extern double Fact[]; //Factorials 1,2,2*3,2*3*4 etc. lookup table.
 extern double Twon[];  //2^n lookup table.
 void InitFactorials();
}
#endif // _Factorial_H_
