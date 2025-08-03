// File: PlaneWaveBS.C  Plane wave basis set.
module;
#include <iosfwd>
#include <vector>

export module qchem.BasisSet.Lattice.PlaneWave;
export import qchem.Lattice;
export import qchem.Types;

import qchem.BasisFunction;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.IBS_Common;
// import qchem.BasisSet.Internal.IEClient;
import qchem.BasisSet.Lattice.Internal.IEClient;

namespace PlaneWave
{

class BasisFunction
    : public Complex_BF
{
public:
    BasisFunction(const RVec3& _k,double _norm) : k(_k), norm(_norm) {};
   
    virtual std::ostream&  Write(std::ostream&   ) const;

    virtual dcmplx  operator()(const RVec3&) const;
    virtual CVec3   Gradient  (const RVec3&) const;


private:
    RVec3 k; //Wave vector, k+G.
    double norm;
};

class IrrepBasisSet
        : public virtual ::IrrepBasisSet,
          public IBS_Common,
          public IrrepIEClient
{
public:
    IrrepBasisSet(IVec3 N, RVec3 k, const std::vector<IVec3>& Gs,double V);
    virtual std::ostream& Write(std::ostream&) const;

};

class BasisSet
    : public BS_Common
{
public:
    BasisSet(const Lattice&, double Emax);
   
};

} //namespace PlaneWave
