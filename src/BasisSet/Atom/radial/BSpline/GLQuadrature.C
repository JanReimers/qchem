// File: GLQuadrature.H Perform Gauss-Legendre quadrature integration over B-Splines.
#include "Imp/BasisSet/Atom/radial/BSpline/GLQuadrature.H"
#include <cassert>

// see gauleg.f 
extern "C"
{
    void gauleg_(const double* rmin, const double* rmax, double* x, double* w, const int* n);
}

GLQuadrature::GLQuadrature(const double& rmin, const double& rmax,int N) 
: xs(N), ws(N) 
{
    gauleg_(&rmin,&rmax,&xs[0],&ws[0],&N); //Numerical recipes.
};


GLCache::GLCache(const bspline::support::Grid<double>& g,size_t N)
: grid(g)
{
    for (size_t i=1;i<g.size();i++)
        itsGLs.push_back(GLQuadrature(g[i-1],g[i],N));
}

double GLCache::Integrate(std::function< double (double)>& f, const sup_t& a, const sup_t& b) const
{
    assert(a.getGrid()==grid);
    assert(b.getGrid()==grid);
    double ret=0;
    sup_t sab=a.calcIntersection(b);
    for (size_t i=0;i<sab.numberOfIntervals();i++)
    {
        size_t ia=sab.absoluteFromRelative(i);
        ret+=itsGLs[ia].Integrate(f);
    }
    return ret;
}

double GLCache::Integrate(std::function< double (double)>& f, const sup_t& a, const sup_t& b, size_t imin, size_t imax) const
{
    assert(imin<imax);
    assert(a.getGrid()==grid);
    assert(b.getGrid()==grid);
    sup_t sab=a.calcIntersection(b);
    if (!sab.containsIntervals()) return 0.0;
    size_t iab_min=sab.absoluteFromRelative(0);
    size_t iab_max=sab.absoluteFromRelative(sab.numberOfIntervals());
    double ret=0;
    for (size_t i=std::max(imin,iab_min);i<std::min(imax,iab_max);i++)
        ret+=itsGLs[i].Integrate(f);
    return ret;
}
