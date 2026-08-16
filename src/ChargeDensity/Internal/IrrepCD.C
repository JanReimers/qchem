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
//! \brief The reciprocal-space trio for a leaf block -- PERIODIC PATH ONLY (V1.7).
//!
//! \c FourierDensity already declared these three pure-virtual, and \c FourierDensityBase<T> already handed
//! it to dcmplx alone -- the mechanism was right there.  The families nevertheless RE-declared the trio in
//! their own bodies for both T, so the finite instantiation had to answer three questions that do not apply
//! to it, and did so with `assert(false)` in an if-constexpr dead branch.  Nine such bodies across the three
//! density families: the largest LSP block in this library.
//!
//! Same CRTP shape as \c IrrepCD_HFPair: the implementation reaches this block's own D and basis without any
//! of it becoming public.
template <class Leaf> class IrrepCD_Fourier : public virtual FourierDensity
{
public:
    //! Metric-free rho-tilde(Delta-m) of THIS block: D contracted against the basis's D-free OVERLAP tensor
    //! Overlap3C(c) (empty kernel) -- the periodic density's native representation.
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const;   // raw rho_DM (0.5(f2)); empty if no raw path
    //! V_H(Delta-m) of THIS block: D contracted against the basis's D-free Coulomb tensor Repulsion3C(c).
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const;
private:
    const Leaf& self() const {return static_cast<const Leaf&>(*this);}
};

//! The reciprocal trio on the periodic path, nothing at all on the finite one.
template <class T, class Leaf> using IrrepFourierBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, IrrepCD_Fourier<Leaf>, NoFourierDensity>;

//! \brief The bra-ket PAIR-PARTNER implementation for a leaf density -- REAL PATH ONLY (V1.6/V1.8).
//!
//! Separated from \c IrrepCD so the complex instantiation never sees it.  Exact exchange is real-only by
//! construction (Vee/Vxc are added by the molecular HF Hamiltonians alone), and the earlier shape -- four
//! virtuals declared on the shared template -- forced \c IrrepCD<dcmplx> to supply four bodies for an
//! operation that does not apply to it.  First they were asserts (silent NO-OPS under -DNDEBUG), then empty
//! stubs; both are the interface failing to segregate.  Now the complex leaf declares nothing at all.
//!
//! CRTP on the leaf so the implementation can use the block's OWN density and basis without any of it
//! becoming public: \c IrrepCD grants friendship, and nothing is exposed through an accessor.
template <class Leaf> class IrrepCD_HFPair : public virtual tHF_Pair_CD<double>
{
public:
    virtual void AccumulateDirectBoth  (rsmat_t& Ji, rsmat_t& Jj, const tHF_Pair_CD<double>& other) const;
    virtual void AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const tHF_Pair_CD<double>& other) const;
    virtual void CompleteDirectPair    (rsmat_t& Ji, rsmat_t& Jj, const rsmat_t& Di,
                                        const BasisSet::Orbital_HF_IBS<double>* bs_i) const;
    virtual void CompleteExchangePair  (rsmat_t& Ki, rsmat_t& Kj, const rsmat_t& Di,
                                        const BasisSet::Orbital_HF_IBS<double>* bs_i) const;
private:
    const Leaf& self() const {return static_cast<const Leaf&>(*this);}
};

//! The pair face on the finite path, nothing at all on the periodic one.
template <class T, class Leaf> using IrrepHF_PairBase =
    std::conditional_t<std::is_same_v<T,double>, IrrepCD_HFPair<Leaf>, NoHF_Pair>;

template <class T> class IrrepCD
    : public virtual tDM_CD<T>
    , public ProjectedDensityBase<T> // AO projection on the finite (double) path; empty on the periodic path
    , public IrrepFourierBase<T,IrrepCD<T>>   // reciprocal trio on the PERIODIC path; NOTHING on the finite one (V1.7)
    , public IrrepHF_PairBase<T,IrrepCD<T>>   // pair partner on the REAL path; NOTHING on the complex one (V1.6)
{
public:
    typedef  mat_t<T>  DenMat;
    typedef hmat_t<T> DenSMat; //Density matrix: HERMITIAN (= symmetric for real T, byte-identical there).

    IrrepCD();
    IrrepCD(const DenSMat&,const tobs_t<T>*, Irrep);

    // The bra-ket pair-partner methods are NOT declared here (V1.6 ISP).  They live in IrrepCD_HFPair,
    // which only the REAL instantiation inherits -- so IrrepCD<dcmplx> declares nothing, needs no vtable
    // slots, and is not forced to define fake bodies for an operation that does not apply to it.
    //! V1.31 whole-system route: this block's basis answers the capability, and the block folds its own
    //! density up to AO.  Real (double) path only -- the dcmplx specializations NA-assert with their siblings.
    virtual const BasisSet::WholeSystemFock_IBS<T>* WholeSystemFock() const;
    virtual void AddAODensity(hmat_t<T>& Dao) const;
    //! AO (auxiliary-basis) projection <rho|c> -- the finite (double) path's ProjectedDensity_AO face; the
    //! periodic (dcmplx) path has no AO face (not cross-cast there), so the dcmplx body is inert.
    virtual double FitGetConstraint() const {return GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const;

    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const;
    virtual rvec_t DM_RhoAtPoints(const rvec3vec_t&, const std::map<Irrep,mat_t<T>>&) const;
    virtual double GetTotalCharge(                      ) const;

    virtual size_t Version() const {return itsVersion;}

    virtual void   ReScale      (double factor              )      ; // No UT coverage
    virtual void   MixIn        (const tMixableDensity<T>&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tMixableDensity<T>&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const rvec3_t&) const;
    //! \f$\nabla\rho\f$ from the density matrix.  Real path only: the \c dcmplx specialization THROWS
    //! (periodic = LDA-only, no complex gradient contraction) rather than returning a silent zero.
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage

    // The reciprocal-space trio is NOT declared here (V1.7 ISP).  It lives in IrrepCD_Fourier, which only
    // the PERIODIC instantiation inherits -- so IrrepCD<double> declares nothing and is not forced to
    // define three "a finite density has no Fourier series" bodies for questions that cannot be put to it.

    virtual std::ostream&       Write(std::ostream&) const;

private:
    friend class IrrepCD_HFPair<IrrepCD<T>>;   // uses this block's own D/basis; nothing is exposed publicly
    friend class IrrepCD_Fourier<IrrepCD<T>>;  // ditto, for the periodic trio
    bool IsZero() const;
    //! The diagonal (self-paired) HF contraction: Jii += J(i,i)·D_i / Kii += K(i,i)·D_i, using THIS block's
    //! own basis on both sides.  Called only from AccumulateDirect/ExchangeBoth's self-pair branch (there is
    //! no bra-ket partner on the diagonal); not part of the tDM_CD interface.
    void AccumulateDirect  (hmat_t<T>& Jii) const;
    void AccumulateExchange(hmat_t<T>& Kii) const;

    DenSMat          itsDensityMatrix;
    const tobs_t<T>* itsBasisSet;
    Spin             itsSpin;
    Irrep        itsIrrep;
    size_t           itsVersion;   //!< TRANSIENT freshness serial (NextDensityVersion); never serialize.
};

} //namespace