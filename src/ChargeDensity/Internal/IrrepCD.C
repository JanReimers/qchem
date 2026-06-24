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
    : public virtual tDM_CD<T>
{
public:
    typedef  mat_t<T>  DenMat;
    typedef hmat_t<T> DenSMat; //Density matrix: HERMITIAN (= symmetric for real T, byte-identical there).

    IrrepCD();
    IrrepCD(const DenSMat&,const tobs_t<T>*, Irrep);

    virtual void AccumulateDirect  (hmat_t<T>& Sab, const ohfbs_t*) const;
    virtual void AccumulateExchange(hmat_t<T>& Sab, const ohfbs_t*) const;
    virtual rvec_t    GetRepulsion3C(const fbs_t*) const;

    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double GetTotalCharge(                      ) const;

    virtual void   ReScale      (double factor              )      ; // No UT coverage
    virtual void   MixIn        (const tDM_CD<T>&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tDM_CD<T>&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const rvec3_t&) const;
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage

    virtual std::ostream&       Write(std::ostream&) const;

private:
    bool IsZero() const;
    
    DenSMat          itsDensityMatrix;
    const tobs_t<T>* itsBasisSet;
    Spin             itsSpin;
    Irrep        itsIrrep;
};

} //namespace