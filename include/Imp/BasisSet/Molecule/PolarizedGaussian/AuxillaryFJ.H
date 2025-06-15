// File: AuxillaryFJ.H   Class for calculating the auxillary function.
#ifndef _AuxillaryFJ_H_
#define _AuxillaryFJ_H_



template <class T> class Vector;

//
//  This class encapsulates the calculation of:
//
//          /1
//  Fj(T) = |  u**(2*j) * exp(-T*u**2)du
//          /0
//
class AuxillaryFJ
{
public:
    void GetFjAt(double T, Vector<double>& Fj) const;
    static const int thejMax;

private:
    static double    theLookUp[121][23];
};

#endif //_AuxillaryFJ_H_
