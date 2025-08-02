// File: BasisSet/Lattice/Internal/IEClient.C Store basis set data needed for integral engines.
module;
#include <vector>
export module qchem.BasisSet.Lattice.Internal.IEClient;
export import qchem.BasisSet.Internal.IEClient;
export import qchem.Types;

export namespace PlaneWave
{

struct IrrepIEClient
    : public virtual ::IrrepIEClient
{
    IrrepIEClient(RVec3 _k,const std::vector<IVec3>& _Gs,double _norm)
        : k(_k), Gs(_Gs), norm(_norm)
        {}

    RVec3 k;
    const std::vector<IVec3>& Gs;
    double norm;

    virtual size_t size() const {return Gs.size();}
};

} //namespace