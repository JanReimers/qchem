// File: PlaneWaveBS.C  Plane wave basis set.
module;
#include <iosfwd>
#include <vector>

export module qchem.BasisSet.Lattice.PlaneWave;
export import qchem.Lattice;
export import qchem.Types;

import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Lattice.Internal.IBS_Evaluator;

namespace PlaneWave
{


class IrrepBasisSet
        : public virtual Complex_IBS
        , public IrrepBasisSet_Common<dcmplx>
        , public IBS_Evaluator
{
public:
    IrrepBasisSet(IVec3 N, RVec3 k, const std::valarray<IVec3>& Gs,double V);
    virtual size_t size() const {return IBS_Evaluator::size();}
    virtual size_t GetNumFunctions() const {return size();}
    virtual const SMatrix<dcmplx>& Overlap() const;
    virtual std::ostream& Write(std::ostream&) const;

};

class BasisSet
    : public BS_Common
{
public:
    BasisSet(const Lattice&, double Emax);
   
};

} //namespace PlaneWave
