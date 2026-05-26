// File: BasisSet1/Atom/Evaluators/BSpline/Internal/GLQuadrature.C Perform Gauss-Legendre quadrature integration over B-Splines.
module;
#include <functional>
#include <bspline/Core.h>
#include <cassert>
#include <iostream>
#include <map>
export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.Types;

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

    template <size_t K> double Integrate(const std::function< double (double)>& w,const bspline::Spline<double,K>& a, const bspline::Spline<double,K>& b) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab,a.getSupport(),b.getSupport());
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

    // used by double Iab_p=gl1.IntegrateIndex(wp,a,b,iab);
    double IntegrateIndex(const std::function< double (double)>& f, size_t i0) const
    {
        assert(i0<itsGLs.size());
        return itsGLs[i0].Integrate(f);
    }
    // used by double Idiag=gl1.IntegrateIndex(wab_diag,a,b,iab); where size_t i1 is a grid index NOT k
    double IntegrateIndex(const std::function< double (double,size_t)>& f, size_t i0) const
    {
        assert(i0<itsGLs.size());
        return itsGLs[i0].Integrate(f);
    }

    // template <class F> double Integrate(const std::function< double (double)>& w,const F& a, const F& b) const
    // {
    //     std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
    //     return Integrate(fwab);
    // }
    template <class F> double Integrate(const std::function< double (double)>& w,const F& a, const F& b, double rmin, double rmax) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab,rmin,rmax);
    }
    size_t RAMsize() const
    {
        size_t ndoubles=grid.size();
        for (auto gl:itsGLs) ndoubles+=gl.RAMsize();
        return ndoubles;
    }
private:
    friend GLCache2D;
    GLCache1D(const GLCache1D&)=delete;
    typedef bspline::Support<double> sup_t;
    double Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b) const;
    double Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b, double rmin, double rmax) const;
    double Integrate(const std::function< double (double)>& f) const
    {
        double ret=0.0;
        for (auto gl:itsGLs)
            ret+=gl.Integrate(f);
        return ret;
    }
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
   
    const bspline::support::Grid<double> grid;
    std::vector<GLQuadrature> itsGLs;
};

export class GLCache2D
{
public:
    GLCache2D(const GLCache1D& gl1,size_t Order);
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
    size_t RAMsize() const;
private:
    GLCache2D(const GLCache2D&)=delete;
    const bspline::support::Grid<double> grid;
    mat_t<GLQuadrature> itsDiagGLs_grid_gl; //Ngrid x Order
    mat_t<GLQuadrature> itsDiagGLs_gl_grid; //Order x Ngrid
};

