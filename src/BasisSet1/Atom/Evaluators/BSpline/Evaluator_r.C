// File: BasisSet1/Atom/Evaluators/BSpline/Evaluator_r.C
module;
#include <bspline/Core.h>
#include <cassert>
#include <functional>
export module qchem.BasisSet.Atom.Evaluators.BSpline.IBS_r;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
import qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet.Atom.Evaluators.Concepts;
import qchem.BasisSet.Internal.Cache4;
// import qchem.Symmetry;
import qchem.Math;

export namespace BasisSet::Atom::Evaluators::BSpline
{
using namespace bspline::integration;
using namespace bspline::operators; 
//
//  This version is for phi(r) = 1/r*sum(Bi(r),i).
// 
template <size_t K> class Evaluator_r : public Internal::EvaluatorCommon<K>, public NR_Angular
{
    using spline_t=Internal::EvaluatorCommon<K>::spline_t;
    using Internal::EvaluatorCommon<K>::splines;
    using Internal::EvaluatorCommon<K>::ns;
    using Internal::EvaluatorCommon<K>::itsGrid;

public: 
    Evaluator_r(size_t Ngrid, double rmin, double rmax, const sym_t& ylm);

    virtual int Getl() const override {return NR_Angular::Getl();}

    double Overlap(size_t i,size_t j) const 
    {
        return BilinearForm{IdentityOperator{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Grad2  (size_t i,size_t j) const 
    {
        static const auto T = -Dx<2>{};
        return BilinearForm{T}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Inv_r1 (size_t i,size_t j) const 
    {
        const spline_t &a=splines[i], &b=splines[j];
        std::function< double (double)> xm1 = [&a,&b](double r)
        {
            assert(r!=0.0);
            return a(r)*b(r)/r;
        };
        return itsGL1D->Integrate(xm1)*FourPi*ns[i]*ns[j];
    } 
    double Inv_r2 (size_t i,size_t j) const 
    {
        const spline_t &a=splines[i], &b=splines[j];
        std::function< double (double)> xm2 = [&a,&b](double r)
        {
            assert(r!=0.0);
            return a(r)*b(r)/(r*r);
        };
        return itsGL1D->Integrate(xm2)*FourPi*ns[i]*ns[j];
    } 
    double Charge (size_t i) const
    {
        return LinearForm{IdentityOperator{}}(splines[i])*ns[i]*FourPi;
    }
    double Norm   (size_t i) const
    {
        return 1.0/sqrt(BilinearForm{IdentityOperator{}}(splines[i],splines[i])*FourPi);
    }

    // Future RKB support.
    // double Grad2  (size_t i,size_t j,const Evaluator_r& b) const 
    // {
    //     static const auto T = -Dx<2>{};
    //     return BilinearForm{T}(splines[i],b.splines[j])*FourPi*ns[i]*b.ns[j]; 
    // } 
    // double Inv_r2 (size_t i,size_t j,const Evaluator_r& _b) const 
    // {
    //     const spline_t &a=splines[i], &b=_b.splines[j];
    //     std::function< double (double)> xm2 = [&a,&b](double r)
    //     {
    //         assert(r!=0.0);
    //         return a(r)*b(r)/(r*r);
    //     };
    //     return itsGL1D->Integrate(xm2)*FourPi*ns[i]*_b.ns[j];
    // } 

    using Evaluators::Evaluator::Norm;
    using Evaluators::Evaluator::size;

    virtual rvec_t     operator() (const rvec3_t&) const override;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const override;

    virtual std::string Name    () const override;
    virtual Cache4*     MakeCache4() const override;
protected:
    rvec_t norms() const; //assumes es,l are already initialized
    std::unique_ptr<GLCache1D> itsGL1D; //We have to hold a pointer, because we don't know grid early enough in the constructor.
};

static_assert(  isOpr_Evaluator<Evaluator_r<6>>);
static_assert(is1E_NR_Evaluator<Evaluator_r<6>>);
static_assert(   isHF_Evaluator<Evaluator_r<6>>);


} //namespace