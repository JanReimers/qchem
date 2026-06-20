// File: AuxillaryFJ.H   Class for calculating the auxillary function.
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.AuxillaryFJ;
import qchem.Types;
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
    void GetFjAt(double T, rvec_t& Fj) const;
    static const int thejMax;

private:
    static double    theLookUp[121][23];
};

