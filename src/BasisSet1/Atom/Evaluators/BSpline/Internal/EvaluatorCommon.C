// File: BasisSet1/Atom/Evaluators/BSpline/Internal/IBS_EvaluatorCommon.C
module;
#include <vector>
#include <bspline/Core.h>

export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
export import qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;

export namespace BasisSet::Atom::Evaluators::BSpline::Internal
{
//
//  This version is for phi(r) = sum(Bi(r),i)
// 
template <size_t K> class EvaluatorCommon : public IBS_Evaluator
{
protected:
    typedef bspline::Spline<double, K> spline_t;
    using rvec11_t=AngularIntegrals::rvec11_t;
public: 
    EvaluatorCommon(size_t Ngrid, double rmin, double rmax, const Irrep_QNs::sym_t& ylm);
    virtual void          Register(Grouper*); //Set up unique spline or exponent indexes.
    virtual std::ostream& Write   (std::ostream&) const;
    virtual size_t        maxSpan () const {return l<=K ? K-l : 0;}  //assume no overlap for indices separated by > maxSpan
    virtual rvec_t        Norm    () const {return ns;}
    static double direct(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::BSpline::RkEngine<K>* cd = dynamic_cast<const ::BSpline::RkEngine<K>*>(c);
        return cd->Coulomb_Rk(la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::BSpline::RkEngine<K>* cd = dynamic_cast<const ::BSpline::RkEngine<K>*>(c);
        return cd->ExchangeRk(la,lc,Ak); // contract over k Rk*Ak, exchange version is more complicated
    }
    virtual std::string RadialID() const;
    virtual std::string RadialType() const;

    const spline_t& operator[](int index) const {return splines[index];}

protected:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);

    using IBS_Evaluator::l;
    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<double>   knots;
    std::vector<spline_t> splines;
    bspline::Grid<double> itsGrid;
};

} //namespace