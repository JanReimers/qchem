// File: BasisSet1/Atom/Evaluators/BSpline/Evaluator.C
module;
#include <bspline/Core.h>
#include <cassert>
export module qchem.BasisSet.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common; 
import qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
import qchem.BasisSet.Atom.Evaluators.Concepts;
import qchem.Symmetry.Spherical;
import qchem.BasisSet.Internal.Cache4;
// import qchem.Symmetry;
import qchem.Math;


export namespace BasisSet::Atom::Evaluators::BSpline
{
using namespace bspline::integration;
using namespace bspline::operators; 
//
//  This version is for phi(r) = sum(Bi(r),i)
// 
template <size_t K> class Evaluator : public Internal::EvaluatorCommon<K>, public NR_Angular
{
    using spline_t=Internal::EvaluatorCommon<K>::spline_t;
    using Internal::EvaluatorCommon<K>::splines;
    using Internal::EvaluatorCommon<K>::ns;
    using Internal::EvaluatorCommon<K>::itsGrid;

public: 
    Evaluator(size_t Ngrid, double rmin, double rmax, const sym_t& ylm);
   
    virtual int Getl() const override {return NR_Angular::Getl();}

    double Overlap(size_t i,size_t j) const 
    {
        return BilinearForm{X<2>{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Grad2  (size_t i,size_t j) const 
    {
        static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
        return BilinearForm{T}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Inv_r1 (size_t i,size_t j) const 
    {
        return BilinearForm{X<1>{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j];
    } 
    double Inv_r2 (size_t i,size_t j) const 
    {
        return BilinearForm{IdentityOperator{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Charge (size_t i) const
    {
        return LinearForm{X<2>{}}(splines[i])*ns[i]*FourPi;
    }
    double Norm   (size_t i) const
    {
        return 1.0/sqrt(BilinearForm{X<2>{}}(splines[i],splines[i])*FourPi);
    }

    // Future RKB support
// double Grad2  (size_t i,size_t j,const Evaluator& b) const 
// {
//     static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
//     return BilinearForm{T}(splines[i],b.splines[j])*FourPi*ns[i]*b.ns[j]; 
// } 
// double Inv_r2 (size_t i,size_t j,const Evaluator& b) const 
// {
//     return BilinearForm{IdentityOperator{}}(splines[i],b.splines[j])*FourPi*ns[i]*b.ns[j]; 
// } 

    using Evaluators::Evaluator::Norm;
    using Evaluators::Evaluator::size;

    virtual rvec_t     operator() (const rvec3_t&) const override;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const override;

    virtual std::string Name      () const override;
    virtual Cache4*     MakeCache4() const override;
protected:
    rvec_t norms() const; //assumes es,l are already initialized
};

static_assert(  isOpr_Evaluator<Evaluator<6>>);
static_assert(is1E_NR_Evaluator<Evaluator<6>>);
static_assert(   isHF_Evaluator<Evaluator<6>>);


} //namespace