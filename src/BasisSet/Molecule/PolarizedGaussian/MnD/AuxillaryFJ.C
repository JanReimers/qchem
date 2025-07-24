// File: AuxillaryFJ.H   Class for calculating the auxillary function.
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.AuxillaryFJ;

import oml;
//
//  This class encapsulates the calculation of:
//
//          /1
//  Fj(T) = |  u**(2*j) * exp(-T*u**2)du
//          /0
//
export class AuxillaryFJ
{
public:
    void GetFjAt(double T, Vector<double>& Fj) const;
    static const int thejMax;

private:
    static double    theLookUp[121][23];
};

