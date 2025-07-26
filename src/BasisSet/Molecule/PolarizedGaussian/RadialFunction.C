// File: RadialFunction.C  Abstract interface for the radial part of a basis function.
module;
#include <vector>
#include <set>
#include <map>
namespace PolarizedGaussian
{

class  Hermite3;
class  RNLM;
}
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;
import Common.UniqueID; 
import qchem.Cluster;
import qchem.ScalarFunction;
import qchem.Streamable;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite1;


// import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;

class  Polarization;

export namespace PolarizedGaussian
{

//
//  Radial functions will be shared by many basis functions, for example 2Px, 2Py and
//  2Pz all have the same radial part, just different polarizations.
//  It is assumed that these radial functions (pointers) will stored in a std::vector 
//  and the SortingIndex is the offset into that list.
//
class RadialFunction
    : public virtual UniqueID
    , public virtual ScalarFunction<double>
    , public virtual Streamable
{
public:
    typedef std::set   <double> sd_t;
    typedef std::vector<double> vd_t;
    virtual       bool      operator==      (const RadialFunction&) const=0; //Ignores L!
    virtual const RVec3&    GetCenter       (                     ) const=0;
    virtual       int       GetL            (                     ) const=0; //N+L+M?
    virtual       double    GetNormalization(const Polarization&  ) const=0;
    virtual       double    GetCharge       (const Polarization&  ) const=0;
    virtual       sd_t      GetExponents    (                     ) const=0; 
    virtual       vd_t      GetCoeff        (                     ) const=0; 

    typedef const RadialFunction rf_t;
    typedef const Polarization   po_t;
    
    virtual double Integrate(qchem::IType2C,rf_t* rb,           po_t& pa, po_t& pb          ,CDCache& cache,const Cluster* cl=0) const=0;
    virtual double Integrate(qchem::IType3C,rf_t* ra, rf_t* rb, po_t& pa, po_t& pb, po_t& pc,CDCache& cache) const=0;
    virtual double Integrate(qchem::IType3C,rf_t* ra,           po_t& pa, po_t& pb, po_t& pc,CDCache& cache, rf_t* rc) const=0;
    
    virtual double Integrate(rf_t* ra,rf_t* rb,rf_t* rc,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache) const=0;
    virtual double Integrate(rf_t* ra,rf_t* rb,         po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache, rf_t* rd) const=0;
    virtual double Integrate(rf_t* ra,                  po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache, rf_t* rc, rf_t* rd) const=0;

    virtual       Hermite3* GetH3 (const RadialFunction&, const RadialFunction&) const=0;
    virtual const Hermite1& GetH1 (                                            ) const=0;

    virtual RadialFunction* Clone      (const RVec3&) const=0;
    virtual RadialFunction* Clone      (            ) const=0;
};

} //namespace PolarizedGaussian

