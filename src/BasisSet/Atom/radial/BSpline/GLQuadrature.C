// File: GLQuadrature.H Perform Gauss-Legendre quadrature integration over B-Splines.
#include "Imp/BasisSet/Atom/radial/BSpline/GLQuadrature.H"
#include <cassert>
#include <iostream>
using std::cout;
using std::endl;

// see gauleg.f 
extern "C"
{
    void gauleg_(const double* rmin, const double* rmax, double* x, double* w, const int* n);
}

GLQuadrature::GLQuadrature(const double& rmin, const double& rmax,int N) 
: its_xmin(rmin), its_xmax(rmax), xs(N), ws(N) 
{
    gauleg_(&rmin,&rmax,&xs[0],&ws[0],&N); //Numerical recipes.
    // cout << "GLQuadrature xmin,xmax,N = " << its_xmin << " " << its_xmax << " " << N << endl;
};


GLCache::GLCache(const bspline::support::Grid<double>& g,size_t N)
: grid(g)
{
    for (size_t i=1;i<grid.size();i++)
        itsGLs.push_back(GLQuadrature(grid[i-1],grid[i],N));
}

double GLCache::Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b) const
{

    assert(a.getGrid()==grid); //Check shared pointers.
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

double GLCache::Integrate(const std::function< double (double)>& f, const sup_t& a, const sup_t& b, double rmin, double rmax) const
{
    if (rmin==rmax) return 0.0;
    if (!(rmin<rmax)) cout << "rmin,rmax=" << rmin << " " << rmax << endl;
    assert(rmin<rmax);
    assert(a.getGrid()==grid);
    assert(b.getGrid()==grid);
    assert(std::isfinite(f(rmin)));
    assert(std::isfinite(f(rmax)));
    sup_t sab=a.calcIntersection(b);
    if (!sab.containsIntervals()) return 0.0; //No support overlap, so interal=0;
    if (rmin<=sab.front() && rmax>=sab.back()) return Integrate(f,a,b); //Integration falls outside support so return full integral.
    auto it_min= std::lower_bound(grid.begin(),grid.end(), rmin); //Get iterator to first grid point *after* rmin
    auto it_max= std::lower_bound(grid.begin(),grid.end(), rmax); //Get iterator to first grid point *after* rmax
    assert(it_min!=grid.end());
    assert(it_max!=grid.end());
    if(it_min!=grid.begin()) it_min--; //Go back one segment
    // cout << grid.front() << " " << *it_min << " " << rmin << " " << rmax << " " << *it_max << " " << grid.back() << endl;
    
    size_t imin= std::distance(grid.begin(), it_min); //These are already absolute.
    size_t imax= std::distance(grid.begin(), it_max);

    size_t iab_min=sab.absoluteFromRelative(0);
    size_t iab_max=sab.absoluteFromRelative(sab.numberOfIntervals());
    double ret=0;
    // cout << iab_min << " " << imin << " " << imax << " " << iab_max << endl;
    for (size_t i=std::max(imin,iab_min);i<std::min(imax,iab_max);i++)
        ret+=itsGLs[i].Integrate(f,rmin,rmax);
    return ret;
}
