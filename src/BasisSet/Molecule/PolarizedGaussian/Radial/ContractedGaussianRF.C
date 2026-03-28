// File: ContractedGaussianRF.H  Radial Contracted Gaussian, no polarization.
module;
#include <vector>
#include <memory>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Radial.ContractedGaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Radial.Common;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite1;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite3;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.Cluster;
import qchem.Types;
import oml;

export namespace PolarizedGaussian
{
//
// Store a list hermite3 pointers.  Contraction coefficients should be absorbed into the scale factors
// of the uncontracted hermite3's before insertion.
//
class ContractedGaussianH3 : public Hermite3
{
public:
    ContractedGaussianH3(const Vector<double>&);
    virtual ~ContractedGaussianH3();

    virtual double operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const;

private:
    friend class ContractedGaussianRF;
    void Insert(Hermite3* h3)
    {
        itsH3s.push_back(std::unique_ptr<Hermite3>(h3));
    }

    typedef std::vector<std::unique_ptr<Hermite3>> h3v_t;
    h3v_t itsH3s;
    const Vector<double>& TheCoeff; //Contraction coefficients.
};

class ContractedGaussianRF
    : public virtual RadialFunction
    , public RadialCommon
{
public:
    ContractedGaussianRF(                                                               );
    ContractedGaussianRF(const Vector<double>& C,std::vector<RadialFunction*>& its_rfs);

    virtual bool            operator==      (const RadialFunction&) const;
    virtual double          GetNormalization(const Polarization&  ) const;
    virtual double          GetCharge       (const Polarization&  ) const;
    virtual sd_t            GetExponents    (                     ) const; 
    virtual vd_t            GetCoeff        (                     ) const;
    
    typedef RadialFunction::rf_t rf_t;
    typedef RadialFunction::po_t po_t;
    virtual double Integrate(       IType  ,rf_t* rb,           po_t& pa, po_t& pb          ,CDCache&,const Cluster* cl=0) const;
    virtual double Integrate(qchem::IType3C,rf_t* ra, rf_t* rb, po_t& pa, po_t& pb, po_t& pc,CDCache&) const;
    virtual double Integrate(qchem::IType3C,rf_t* ra,           po_t& pa, po_t& pb, po_t& pc,CDCache&, rf_t* rc) const;

    virtual double Integrate(rf_t* ra,rf_t* rb,rf_t* rc,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&) const;
    virtual double Integrate(rf_t* ra,rf_t* rb,         po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&, rf_t* rd) const;
    virtual double Integrate(rf_t* ra,                  po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&, rf_t* rc, rf_t* rd) const;

    
    virtual Hermite3*       GetH3 (const RadialFunction&, const RadialFunction&) const;

    virtual std::ostream&        Write(std::ostream&    ) const;
    virtual RadialFunction* Clone(            ) const;
    virtual RadialFunction* Clone(const rvec3_t&) const;

    virtual double          operator()(const rvec3_t&) const;
    virtual rvec3_t           Gradient  (const rvec3_t&) const;

private:
    virtual Hermite1* MakeH1() const;

    void Check() const;

    typedef std::vector<std::unique_ptr<RadialFunction>> rfv_t;

    Vector<double> cs;
    Vector<double> unormalized_cs;
    rfv_t          gs;
};

} //namespace PolarizedGaussian

