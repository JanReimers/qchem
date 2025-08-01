// File: GaussianH3.C  Class for managing 3 function Hermite coefficients.
module;
#include <iosfwd>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Radial.GaussianH3;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite3;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.Types;
import oml; 

#define LMAX 3
export namespace PolarizedGaussian
{
//
//  Encapsulate calculation and storage of Hermite expansion coefficients
//    _=
//   nnn
//  d
//   0
//
//  No attempt at optimizing storage space is made here as these data structures
//  are intended for "on the fly" calculation.
//  The Scale parameter is applied to every data element, if scale=Eabc*(Pi/AlphaQ)^3/2 the
//  op() returns the overlap integral directly.
//


class GaussianH3 : public Hermite3
{
public:
    GaussianH3();
    GaussianH3(double AlphaP, const RVec3& PA, const RVec3& PB, const RVec3& PC, int LA, int LB, int LC, double Scale=1.0);

    virtual double operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const;


private:
    typedef double Array4D[3*LMAX+1][LMAX+1][LMAX+1][LMAX+1];
    friend std::ostream& operator<<(std::ostream&,const GaussianH3&);

    double GetAny(const Array4D, int N, int na, int nb, int nc) const;
    double Getd(int N, int na, int nb, int nc) const
    {
        return GetAny(d,N,na,nb,nc);
    }
    double Gete(int N, int na, int nb, int nc) const
    {
        return GetAny(e,N,na,nb,nc);
    }
    double Getf(int N, int na, int nb, int nc) const
    {
        return GetAny(f,N,na,nb,nc);
    }

    double   itsa12s[3*LMAX+1];
    Array4D  d,e,f;
    int      itsLA, itsLB, itsLC;
    double   itsScale;
};

} //namespace PolarizedGaussian

