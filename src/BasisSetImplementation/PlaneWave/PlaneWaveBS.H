// File: PlaneWaveBS.H  Polarized Gaussian basis set, for MO calculations.
#ifndef _PlaneWaveBS_H_
#define _PlaneWaveBS_H_



#include "BasisSetImplementation/BasisSetImplementation.H"
#include "BasisSetImplementation/TBasisSetImplementation.H"
#include "Cluster/UnitCell.H"

class Reader;
class Cluster;
class Lattice;

class PlaneWaveBS
    : public virtual TIrrepBasisSet<std::complex<double> >
    , public         BasisSetImplementation
    , public         TBasisSetImplementation<std::complex<double> >
{
public:
    PlaneWaveBS();
    PlaneWaveBS(IntegralDataBase<std::complex<double> >*, const Cluster*, double Gmax, Mesh* m=0);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual IrrepBasisSet* Clone(            ) const;
    virtual IrrepBasisSet* Clone(const RVec3&) const;

private:
    void MakeBasisFunctions();

    friend class PlaneWaveIE;

    UnitCell             itsRLCell; //Unit cell for the reciprocal lattice (scaled by 2*Pi).
};

#endif //_PlaneWaveBS_H_
