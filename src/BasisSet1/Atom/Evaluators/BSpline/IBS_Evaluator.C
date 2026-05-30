// File: BasisSet1/Atom/Evaluators/BSpline/IBS_Evaluator.C
module;
#include <bspline/Core.h>
#include <cassert>
#include <cmath>
export module qchem.BasisSet.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
// import qchem.BasisSet.Atom.Evaluators.IBS;
// import qchem.BasisSet.Atom.Evaluators.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.Concepts;
import qchem.Symmetry.Yl;
import Common.Constants;
import qchem.BasisSet.Internal.Cache4;
import Common.IntPow;

// required by BSpline_Cache4
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;

export namespace BasisSet::Atom::Evaluators::BSpline
{
    using namespace bspline::integration;
    using namespace bspline::operators; 
//
//  This version is for phi(r) = sum(Bi(r),i)
// 
template <size_t K> class BSpline_IBS_Evaluator : public Internal::EvaluatorCommon<K>
{
    // typedef bspline::Spline<double, K> spline_t;
    using spline_t=Internal::EvaluatorCommon<K>::spline_t;
    using Internal::EvaluatorCommon<K>::splines;
    using Internal::EvaluatorCommon<K>::ns;
    using Internal::EvaluatorCommon<K>::itsGrid;
    using Internal::EvaluatorCommon<K>::l;

public: 
    BSpline_IBS_Evaluator(size_t Ngrid, double rmin, double rmax, const Irrep_QNs::sym_t& ylm);

    double Overlap(size_t i,size_t j) const //no l dependence
    {
        return BilinearForm{X<2>{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Grad2  (size_t i,size_t j) const //no l dependence
    {
        static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
        return BilinearForm{T}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Grad2  (size_t i,size_t j,const BSpline_IBS_Evaluator& b) const //no l dependence
    {
        static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
        return BilinearForm{T}(splines[i],b.splines[j])*FourPi*ns[i]*b.ns[j]; 
    } 
    double Inv_r1 (size_t i,size_t j) const //no l dependence
    {
        return BilinearForm{X<1>{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j];
    } 
    double Inv_r2 (size_t i,size_t j) const //no l dependence
    {
        return BilinearForm{IdentityOperator{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Inv_r2 (size_t i,size_t j,const BSpline_IBS_Evaluator& b) const //no l dependence
    {
        return BilinearForm{IdentityOperator{}}(splines[i],b.splines[j])*FourPi*ns[i]*b.ns[j]; 
    } 
    double Charge (size_t i) const
    {
        return LinearForm{X<2>{}}(splines[i])*ns[i]*FourPi;
    }
    double Norm   (size_t i) const
    {
        return 1.0/sqrt(BilinearForm{X<2>{}}(splines[i],splines[i])*FourPi);
    }

    using IBS_Evaluator::Norm;

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name    () const;
    virtual Cache4*    MakeCache4() const;
protected:
    rvec_t norms() const; //assumes es,l are already initialized
};

static_assert(isGeneric_Evaluator<BSpline_IBS_Evaluator<6>>);
static_assert(is1E_Evaluator     <BSpline_IBS_Evaluator<6>>);
// static_assert(isFit_Evaluator    <BSpline_IBS_Evaluator<6>>);
// static_assert(isDFT_Evaluator    <BSpline_IBS_Evaluator<6>>);
static_assert(isRKBL_Evaluator   <BSpline_IBS_Evaluator<6>>);
static_assert(isHF_Evaluator     <BSpline_IBS_Evaluator<6>>);

std::ostream& operator<<(std::ostream& os, const bspline::Support<double>& sup)
{
    return os << "[" << sup.front() << "," << sup.back() << "]";
}

template <size_t K> class BSpline_Cache4 : public  Cache4
{
public:
    BSpline_Cache4(const bspline::Grid<double>& grid);
    ~BSpline_Cache4() {delete itsRkCache;}
    virtual void   Register(Cache4_Client * eval);
    virtual Rk*    Create  (size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual size_t RAMsize () const;

private:
    using func_t=::BSpline::RkEngine<K>::func_t;
    func_t wp,wm; //Weight functions for Slater integrals.
    size_t itsMaxl;
    GLCache1D   itsGL1D;
    GLCache2D   itsGL2D;

    SplineGrouper<K> grouper;
    ::BSpline::RkCache<K>* itsRkCache;
};

} //namespace