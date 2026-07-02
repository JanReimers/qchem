// File: BasisSet/Atom/Evaluators/BSpline/Internal/EvaluatorCommon.C
module;
#include <vector>
#include <memory>
#include <utility>
#include <bspline/Core.h>

export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
export import qchem.BasisSet.Atom.Evaluators;
export import qchem.Symmetry;
import qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.Internal.Grouper;
// required by Cache4
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;

export namespace qchem::BasisSet::Atom::Evaluators::BSpline::Internal
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
    static double direct(const Cacheable4* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::qchem::BSpline::RkEngine<K>* cd = dynamic_cast<const ::qchem::BSpline::RkEngine<K>*>(c);
        return cd->DirectRk  (la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable4* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::qchem::BSpline::RkEngine<K>* cd = dynamic_cast<const ::qchem::BSpline::RkEngine<K>*>(c);
        return cd->ExchangeRk(la,lc,Ak); // contract over k Rk*Ak, exchange version is more complicated
    }
    virtual std::ostream& Write(std::ostream&) const override;
    virtual std::string   RadialID  () const override;
    virtual std::string   RadialType() const override;

    const spline_t& operator[](int index) const {return splines[index];}
    const bspline::Grid<double>& GetGrid() const {return itsGrid;}

protected:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax, int l);

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<double>   knots;
    std::vector<spline_t> splines;
    bspline::Grid<double> itsGrid;
    rvec_t ns;
    std::vector<size_t> es_indices; //Unique spline index
};

// Grid-specific quadrature + Rk-moment tables.  A single "BSpline<K>" Cache4 now serves EVERY grid of
// order K (RadialType()==Name(), like SG/SL), so the grid-dependent Gauss-Legendre caches live here --
// one bundle per distinct grid -- instead of forcing a whole separate Cache4 (and grid string in the key)
// per grid.
template <size_t K> struct GridData
{
    GridData(const bspline::Grid<double>& g, size_t Kp) : grid(g), gl1(g,K+Kp), gl2(g,2*K+Kp,K+3) {}
    bspline::Grid<double> grid;
    GLCache1D gl1;
    GLCache2D gl2;
    mutable std::unique_ptr<::qchem::BSpline::RkCache<K>> rkcache; //built lazily on first Create (see below)
};

template <size_t K> class Cache4 : public  ::qchem::Cache4
{
    using func_t=::qchem::BSpline::RkEngine<K>::func_t;
public:
    Cache4(const func_t& wp, const func_t& wm, size_t Kp);
    virtual void   Register(Cache4_Client * eval);
    virtual Rk*    Create  (size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual size_t RAMsize () const;

private:
    GridData<K>&       ensureGrid(const bspline::Grid<double>& g);       //find-or-build the per-grid bundle
    const GridData<K>& gridFor   (const bspline::Grid<double>& g) const; //Create() only asks for built grids

    func_t wp,wm; //Weight functions for Slater integrals.
    size_t Kp;
    size_t itsMaxl;
    SplineGrouper<K> grouper; //lossless across grids -- separates splines from different knot vectors
    std::vector<std::unique_ptr<GridData<K>>> itsGrids; //one per distinct grid, keyed by GridData::grid
};

} //namespace