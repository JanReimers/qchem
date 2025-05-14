// File: GLQuadrature.H Perform Gauss-Legendre quadrature integration over B-Splines.
#ifndef _GLQuadrature_H_
#define _GLQuadrature_H_

#include <valarray>
#include <functional>
#include <bspline/Core.h>

class GLQuadrature
{
public:
    GLQuadrature(const double& rmin, const double& rmax,int N);

    double Integrate(std::function< double (double)>& f) const
    {
        double ret=0.0;
        for (size_t i=0;i<xs.size();i++)
            ret+=ws[i]*f(xs[i]);
        return ret;
    }
private:
    std::valarray<double> xs,ws;
};

class GLCache
{
public:
    GLCache(const bspline::support::Grid<double>& g,size_t N);
    
    template <size_t K> double Integrate(std::function< double (double)>& w,const bspline::Spline<double,K>& a, const bspline::Spline<double,K>& b) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab,a.getSupport(),b.getSupport());
    }
    template <size_t K> double Integrate(std::function< double (double)>& w,const bspline::Spline<double,K>& a, const bspline::Spline<double,K>& b, size_t imin, size_t imax) const
    {
        std::function< double (double)> fwab = [w,a,b](double x){return w(x)*a(x)*b(x);};
        return Integrate(fwab,a.getSupport(),b.getSupport(),imin,imax);
    }

private:
    typedef bspline::Support<double> sup_t;
    double Integrate(std::function< double (double)>& f, const sup_t& a, const sup_t& b) const;
    double Integrate(std::function< double (double)>& f, const sup_t& a, const sup_t& b, size_t imin, size_t imax) const;
    const bspline::support::Grid<double>& grid;
    std::vector<GLQuadrature> itsGLs;
};


#endif //_GLQuadrature_H_