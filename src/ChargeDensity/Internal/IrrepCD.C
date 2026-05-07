// File: ExactIrrepCD.H  Exact charged density for ONE irreducable representation basis set.
module;
#include <iosfwd>
export module qchem.ChargeDensity.Imp.IrrepCD;

export import qchem.ChargeDensity;
export import qchem.Symmetry.Irrep;
import qchem.ChargeDensity.Types;

export namespace qchem::ChargeDensity
{

//------------------------------------------------------------------------------------
//
//  This maintains the exact charge density represented by the density matrix
//  of one irreducable representation.  The full charge density will in general
//  be a summation of these guys.
//
template <class T> class IrrepCD
    : public virtual DM_CD
{
public:
    typedef  mat_t<T>  DenMat;
    typedef smat_t<T> DenSMat; //Type for the density matrix.
    
    IrrepCD();
    IrrepCD(const DenSMat&,const tobs_t<T>*, Irrep_QNs);

    virtual void AccumulateDirect  (rsmat_t& Sab, const ohfbs_t*) const;
    virtual void AccumulateExchange(rsmat_t& Sab, const ohfbs_t*) const;
    virtual rvec_t    GetRepulsion3C(const Fit_IBS*) const;

    virtual double DM_Contract(const Static_CC*) const;
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const;
    virtual double GetTotalCharge(                      ) const;

    virtual void   ReScale      (double factor         )      ; // No UT coverage
    virtual void   MixIn        (const DM_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const rvec3_t&) const;
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage

    virtual std::ostream&       Write(std::ostream&) const;

private:
    bool IsZero() const;
    
    DenSMat          itsDensityMatrix;
    const tobs_t<T>* itsBasisSet;
    Spin             itsSpin;
    Irrep_QNs        itsIrrep;
};

} //namespace