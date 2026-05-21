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
    GLQuadrature() {}; //map needs a default constructor.
    GLQuadrature(const double& rmin, const double& rmax,int N);

    double Integrate(const std::function< double (double)>& f) const
    {
        double ret=0.0;
        for (size_t i=0;i<xs.size();i++)
            ret+=ws[i]*f(xs[i]);
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

    template <size_t K1, size_t K2> double Integrate(const std::function< double (double)>& w,const bspline::Spline<double,K1>& a, const bspline::Spline<double,K2>& b) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab,a.getSupport(),b.getSupport());
    }
    template <size_t K> double Integrate(const std::function< double (double)>& w,const bspline::Spline<double,K>& a, const bspline::Spline<double,K>& b) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab,a.getSupport(),b.getSupport());
    }
    template <size_t K> double Integrate(const std::function< double (double)>& w,const bspline::Spline<double,K>& a, const bspline::Spline<double,K>& b, double rmin, double rmax) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab,a.getSupport(),b.getSupport(),rmin,rmax);
    }
    template <size_t K> double IntegrateIndex(const std::function< double (double)>& w,const bspline::Spline<double,K>& a, const bspline::Spline<double,K>& b, size_t i0, size_t i1) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return IntegrateIndex(fwab,i0,i1);
    }
    template <size_t K> double IntegrateIndex(const std::function< double (double)>& w,const bspline::Spline<double,K>& a, const bspline::Spline<double,K>& b, size_t i0) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        assert(i0<itsGLs.size());
        return itsGLs[i0].Integrate(fwab);
    }

    template <class F> double Integrate(const std::function< double (double)>& w,const F& a, const F& b) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab);
    }
    template <class F> double Integrate(const std::function< double (double)>& w,const F& a, const F& b, double rmin, double rmax) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab,rmin,rmax);
    }

private:
    friend GLCache2D;
    GLCache1D(const GLCache1D&)=delete;
    typedef bspline::Support<double> sup_t;
    double Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b) const;
    double Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b, double rmin, double rmax) const;
    double IntegrateIndex(const std::function< double (double)>& f,  size_t i0, size_t i1) const
    {
        if (i0>=i1) return 0.0;
        assert(i0<itsGLs.size());
        assert(i1<=itsGLs.size());
        double ret=0.0;
        for (size_t i=i0;i<i1;i++) ret+=itsGLs[i].Integrate(f);
        return ret;
    }
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
    const GLQuadrature& find(double rmin, double rmax) const;
private:
    GLCache2D(const GLCache2D&)=delete;
    const bspline::support::Grid<double> grid;
    std::map<double,std::map<double,GLQuadrature>> itsDiagGLs;

};

