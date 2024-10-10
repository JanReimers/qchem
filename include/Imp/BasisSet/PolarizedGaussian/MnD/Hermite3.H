// File: Hermite3.H  Interface class for managing 3 function Hermite coefficients
#ifndef _Hermite3_H_
#define _Hermite3_H_



//
//  Encapsulate calculation and storage of Hermite expansion coefficients
//    _=
//   nnn
//  d
//   0
//
//  This is an interface base class for primative and contracted Hermite3 blocks.
//

class Polarization;

class Hermite3
{
public:
    virtual ~Hermite3() {};
    virtual double operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const =0;
};

#endif //_Hermite3_H_
