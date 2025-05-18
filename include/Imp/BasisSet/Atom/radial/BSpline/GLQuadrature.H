// File: GLQuadrature.H Perform Gauss-Legendre quadrature integration over B-Splines.
#ifndef _GLQuadrature_H_
#define _GLQuadrature_H_

#include <valarray>
#include <functional>
#include <bspline/Core.h>
#include <cassert>
#include <iostream>

class GLQuadrature
{
public:
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
    double its_xmin, its_xmax;
    std::valarray<double> xs,ws;
};

class GLCache
{
public:
    GLCache(const bspline::support::Grid<double>& g,size_t N);
    
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

// private:
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


#endif //_GLQuadrature_H_