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
        // std::cout << "Xmin,max(xs) = " << xmin << " " << xs[xs.size()-1] << std::endl;
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

private:
    typedef bspline::Support<double> sup_t;
    double Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b) const;
    double Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b, double rmin, double rmax) const;
    const bspline::support::Grid<double>& grid;
    std::vector<GLQuadrature> itsGLs;
};


#endif //_GLQuadrature_H_