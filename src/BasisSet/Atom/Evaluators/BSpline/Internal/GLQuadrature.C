// File: BasisSet/Atom/Evaluators/BSpline/Internal/GLQuadrature.C Perform Gauss-Legendre quadrature integration over B-Splines.
module;
#include <functional>
#include <bspline/Core.h>
#include <cassert>
#include <iostream>
#include <map>
export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.Types;

namespace qchem {

export class GLCache1D;
export class GLCache2D;

export class GLQuadrature
{
public:
    GLQuadrature() : xs(0), ws(0) {}; //map needs a default constructor.
    GLQuadrature(const double& rmin, const double& rmax,int N);

    double Integrate(const std::function< double (double)>& f) const
    {
        double ret=0.0;
        for (size_t i=0;i<xs.size();i++)
            ret+=ws[i]*f(xs[i]);
        return ret;
    }
    double Integrate(const std::function< double (double, size_t)>& f) const
    {
        double ret=0.0;
        for (size_t i=0;i<xs.size();i++)
            ret+=ws[i]*f(xs[i],i);
        return ret;
    }
    //! Weighted sum of a PRE-SAMPLED product integrand: sum_i weight(i) * ( w(node(i),k) * u[i] * v[i] ).
    //! The caller samples the two factors u=f(node),v=g(node) ONCE (e.g. two B-splines in the Rk 2D
    //! quadrature) and reuses them across all k; this reproduces Integrate(f) with f=w*u*v term-for-term
    //! (same order and weight*(w*u*v) associativity) but never re-evaluates the splines.  Templated on the
    //! sample containers so blaze row/column views pass through with no std::function indirection.
    template <class W, class U, class V>
    double Integrate(const W& w, size_t k, const U& u, const V& v) const
    {
        double ret=0.0;
        for (size_t i=0;i<xs.size();i++)
            ret+=ws[i]*(w(xs[i],k)*u[i]*v[i]);
        return ret;
    }
    double Integrate(const std::function< double (double)>& f, double xmin, double xmax) const
    {
        assert(xmax>=its_xmin); //Make sure caller did thier homework and checked the ranges.
        assert(xmin<=its_xmax);
        double ret=0.0;
        for (size_t i=0;i<xs.size();i++)
            if (xs[i]>=xmin && xs[i]<=xmax)
                ret+=ws[i]*f(xs[i]);
        return ret;
    }
    size_t RAMsize() const
    {
        return 2+2*xs.size();
    }
    //! Node access so a caller can pre-sample a fixed integrand (e.g. a B-spline) at the quadrature points
    //! ONCE and then reuse it across many weight functions via the Integrate(w,k,u,v) overload below --
    //! avoiding the per-call spline re-evaluation (and its knot-span binary search).
    size_t size()           const {return xs.size();}
    double node(size_t i)   const {return xs[i];}
private:
    friend class GLCache1D;
    friend class GLCache2D;
    double its_xmin, its_xmax;
    rvec_t xs,ws;
};
export class GLCache1D
{
public:
    GLCache1D(const bspline::support::Grid<double>& g,size_t Order);
    template <size_t K> using sp_t=bspline::Spline<double,K>;
    template <size_t K> double Integrate(const std::function< double (double)>& w, const sp_t<K>& a, const sp_t<K>& b ) const
    {
        std::function< double (double)> wab=[&w,&a,&b](double x) {return w(x)*a(x)*b(x);};
        return Integrate(wab,a.getSupport(),b.getSupport());
    }
    template <size_t K> double Integrate(const std::function< double (double)>& w, const sp_t<K>& a, const sp_t<K>& b , double rmin, double rmax) const
    {
        std::function< double (double)> wab=[&w,&a,&b](double x) {return w(x)*a(x)*b(x);};
        return Integrate(wab,a.getSupport(),b.getSupport(),rmin,rmax);
    }
    template <size_t K> double Integrate(const std::function< double (double,size_t)>& w, const sp_t<K>& a, const sp_t<K>& b, size_t k ) const
    {
        std::function< double (double)> wab=[k,&w,&a,&b](double x) {return w(x,k)*a(x)*b(x);};
        return Integrate(wab,a.getSupport(),b.getSupport());
    }
    ///used by double Icd_p=gl1.IntegrateIndex(wpcd,scd.getStartIndex(),iab);
    double IntegrateIndex(const std::function< double (double)>& f,  size_t i0, size_t i1) const
    {
        if (i0>=i1) return 0.0;
        assert(i0<itsGLs.size());
        assert(i1<=itsGLs.size());
        double ret=0.0;
        for (size_t i=i0;i<i1;i++) ret+=itsGLs[i].Integrate(f);
        return ret;
    }
    // Used in some unit tests.
    double Integrate(const std::function< double (double)>& f, double rmin, double rmax) const
    {
        double ret=0.0;
        // std::cout << grid.size() << " " << itsGLs.size() << std::endl;
        for (size_t i=1;i<grid.size();i++)
        {
            // std::cout << "rmin,grid[i-1],grid[i],rmax=" << rmin << "   " <<  grid[i-1] << "   " << grid[i] << "   " << rmax << std::endl;
            if (rmin<=grid[i] && rmax>=grid[i-1])
                ret+=itsGLs[i-1].Integrate(f,rmin,rmax);

        }
        return ret;
    }
    double Integrate(const std::function< double (double)>& f) const
    {
        double ret=0.0;
        for (auto gl:itsGLs)
            ret+=gl.Integrate(f);
        return ret;
    }
    const GLQuadrature& operator[](size_t index) const 
    {
        assert(index<itsGLs.size());
        return itsGLs[index];
    }
    size_t RAMsize() const
    {
        size_t ndoubles=grid.size();
        for (auto gl:itsGLs) ndoubles+=gl.RAMsize();
        return ndoubles;
    }
private:
    typedef bspline::Support<double> sup_t;
    double Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b) const;
    double Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b, double rmin, double rmax) const;
    
    friend GLCache2D;
    GLCache1D(const GLCache1D&)=delete;
    const bspline::support::Grid<double> grid;
    std::vector<GLQuadrature> itsGLs;
};

export class GLCache2D
{
public:
    // GLCache2D(const GLCache1D& gl1,size_t Order);
    GLCache2D(const bspline::support::Grid<double>& g,size_t Order1, size_t Order2);
    const GLQuadrature& find_grid_gl(size_t igrid,size_t igl) const
    {
        assert(igrid>=0);
        assert(igl>=0);
        assert(igrid<itsDiagGLs_grid_gl.rows());
        assert(igl<itsDiagGLs_grid_gl.columns());
        return itsDiagGLs_grid_gl(igrid,igl);
    }
    const GLQuadrature& find_gl_grid(size_t igl,size_t igrid) const
    {
        assert(igrid>=0);
        assert(igl>=0);
        assert(igl<itsDiagGLs_gl_grid.rows());
        assert(igrid<itsDiagGLs_gl_grid.columns());
        return itsDiagGLs_gl_grid(igl,igrid);
    }
    const GLCache1D& GL1d() const {return itsGl1D;}
    size_t RAMsize() const;
private:
    GLCache2D(const GLCache2D&)=delete;
    const bspline::support::Grid<double> grid;
    GLCache1D itsGl1D;
    mat_t<GLQuadrature> itsDiagGLs_grid_gl; //Ngrid x Order
    mat_t<GLQuadrature> itsDiagGLs_gl_grid; //Order x Ngrid
};


} // namespace qchem