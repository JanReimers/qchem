// File: BasisSet1/Atom/Evaluators/BSpline/Internal/IBS_EvaluatorCommon.C
module;
#include <vector>
#include <bspline/Core.h>

export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
export import qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;

// required by BSpline_Cache4
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;

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

template <size_t K> class BSpline_Cache4 : public  ::Cache4
{
    using func_t=::BSpline::RkEngine<K>::func_t;
public:
    BSpline_Cache4(const bspline::Grid<double>& grid,const func_t& wp, const func_t& wm);
    ~BSpline_Cache4() {delete itsRkCache;}
    virtual void   Register(Cache4_Client * eval);
    virtual Rk*    Create  (size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual size_t RAMsize () const;

private:
    func_t wp,wm; //Weight functions for Slater integrals.
    size_t itsMaxl;
    GLCache1D   itsGL1D;
    GLCache2D   itsGL2D;

    SplineGrouper<K> grouper;
    ::BSpline::RkCache<K>* itsRkCache;
};

} //namespace