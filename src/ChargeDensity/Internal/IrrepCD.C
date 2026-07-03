// File: ExactIrrepCD.H  Exact charged density for ONE irreducable representation basis set.
module;
#include <iosfwd>
#include <cstddef>
#include <map>
#include <string>
export module qchem.ChargeDensity.Imp.IrrepCD;

export import qchem.ChargeDensity;
export import qchem.Symmetry.Irrep;
export import qchem.ChargeDensity.FourierDensity;   // G-space rho-tilde (plane-wave path)
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
    , public ProjectedDensityBase<T> // AO projection on the finite (double) path; empty on the periodic path
    , public FourierDensityBase<T>   // FourierDensity on the periodic (dcmplx) path; empty on the finite path
{
public:
    typedef  mat_t<T>  DenMat;
    typedef hmat_t<T> DenSMat; //Density matrix: HERMITIAN (= symmetric for real T, byte-identical there).

    IrrepCD();
    IrrepCD(const DenSMat&,const tobs_t<T>*, Irrep);

    virtual void AccumulateDirect  (hmat_t<T>& Sab, const ohfbs_t*) const;
    virtual void AccumulateExchange(hmat_t<T>& Sab, const ohfbs_t*) const;
    //! Bra-ket pair partner (doc/ERI4Rework.md §4): scatter the canonical block J(this,other) into BOTH
    //! Fock blocks Ji,Jj.  \a other is a sibling IrrepCD (same-class cast, as MixIn/GetChangeFrom do).
    virtual void AccumulateDirectBoth  (hmat_t<T>& Ji, hmat_t<T>& Jj, const tDM_CD<T>& other) const;
    virtual void AccumulateExchangeBoth(hmat_t<T>& Ki, hmat_t<T>& Kj, const tDM_CD<T>& other) const;
    //! AO (auxiliary-basis) projection <rho|c> -- the finite (double) path's ProjectedDensity_AO face; the
    //! periodic (dcmplx) path has no AO face (not cross-cast there), so the dcmplx body is inert.
    virtual double FitGetConstraint() const {return GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::FIT_CD_ABS*) const;

    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const;
    virtual double GetTotalCharge(                      ) const;

    virtual size_t Version() const {return itsVersion;}

    virtual void   ReScale      (double factor              )      ; // No UT coverage
    virtual void   MixIn        (const tDM_CD<T>&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tDM_CD<T>&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const rvec3_t&) const;
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage

    //! Reciprocal-space coefficients rho-tilde(Delta-m) of THIS block (= basis->MakeFourierDensity(D)).
    //! The periodic density's native representation; a finite density has none (real path NA-asserts).
    virtual FourierMap GetFourierDensity() const;

    virtual std::ostream&       Write(std::ostream&) const;

private:
    bool IsZero() const;
    
    DenSMat          itsDensityMatrix;
    const tobs_t<T>* itsBasisSet;
    Spin             itsSpin;
    Irrep        itsIrrep;
    size_t           itsVersion;   //!< TRANSIENT freshness serial (NextDensityVersion); never serialize.
};

} //namespace