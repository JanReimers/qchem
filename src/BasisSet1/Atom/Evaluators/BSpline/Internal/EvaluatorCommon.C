// File: BasisSet1/Atom/Evaluators/BSpline/Internal/EvaluatorCommon.C
module;
#include <vector>
#include <bspline/Core.h>

export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
export import qchem.BasisSet.Atom.Evaluators;
export import qchem.Symmetry;
import qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
// required by Cache4
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;

export namespace BasisSet::Atom::Evaluators::BSpline::Internal
{
std::ostream& operator<<(std::ostream& os, const bspline::Support<double>& sup)
{
    return os << "[" << sup.front() << "," << sup.back() << "]";
}
std::ostream& operator<<(std::ostream& os, const bspline::Grid<double>& grid)
{
    os << "{";
    for (auto g:grid) os << g << ",";

    return os << "}";
}


template <size_t K> class EvaluatorCommon : public virtual Evaluators::Evaluator , public virtual Evaluators::HF_Evaluator
{
protected:
    typedef bspline::Spline<double, K> spline_t;
    using rvec11_t=rvec11_t;
public: 
    EvaluatorCommon(size_t Ngrid, double rmin, double rmax, const sym_t& ylm);
    virtual void   Register(Grouper*) override; //Set up unique spline or exponent indexes.
    virtual size_t maxSpan () const override {return Getl()<=K ? K-Getl() : 0;}  //assume no overlap for indices separated by > maxSpan
    virtual size_t size    () const override { return ns.size(); }
    virtual size_t es_index(size_t i     ) const {return es_indices[i];}
    virtual rvec_t Norm    () const override { return ns; }
    static double direct(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::BSpline::RkEngine<K>* cd = dynamic_cast<const ::BSpline::RkEngine<K>*>(c);
        return cd->DirectRk  (la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::BSpline::RkEngine<K>* cd = dynamic_cast<const ::BSpline::RkEngine<K>*>(c);
        return cd->ExchangeRk(la,lc,Ak); // contract over k Rk*Ak, exchange version is more complicated
    }
    virtual std::ostream& Write(std::ostream&) const override;
    virtual std::string   RadialID  () const override;
    virtual std::string   RadialType() const override;

    const spline_t& operator[](int index) const {return splines[index];}

protected:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax, int l);

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<double>   knots;
    std::vector<spline_t> splines;
    bspline::Grid<double> itsGrid;
    rvec_t ns;
    std::vector<size_t> es_indices; //Unique spline index
};

template <size_t K> class Cache4 : public  ::Cache4
{
    using func_t=::BSpline::RkEngine<K>::func_t;
public:
    Cache4(const bspline::Grid<double>& grid,const func_t& wp, const func_t& wm, size_t Kp);
    ~Cache4() {delete itsRkCache;}
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