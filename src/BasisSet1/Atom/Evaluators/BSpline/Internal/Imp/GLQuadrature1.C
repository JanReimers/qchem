// File: GLQuadrature.H Perform Gauss-Legendre quadrature integration over B-Splines.
module;
#include <cassert>
#include <iostream>
#include <functional>
#include <map>
#include <bspline/Core.h>
#include <blaze/Math.h>

module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
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
}

GLCache1D::GLCache1D(const bspline::support::Grid<double>& g,size_t Order)
: grid(g)
{
    for (size_t i=1;i<grid.size();i++)
        itsGLs.push_back(GLQuadrature(grid[i-1],grid[i],Order));
}

GLCache2D::GLCache2D(const GLCache1D& gl1,size_t Order)
: grid(gl1.grid), itsDiagGLs_grid_gl(grid.size(),Order), itsDiagGLs_gl_grid(Order,grid.size())
{
    for (size_t i=1;i<grid.size();i++)
    {
        double rmin=grid[i-1], rmax=grid[i];
        size_t j=0;
        for (double r:gl1.itsGLs[i-1].xs)
        {
            itsDiagGLs_grid_gl(i-1,j)=GLQuadrature(rmin,r,Order);
            itsDiagGLs_gl_grid(j  ,i)=GLQuadrature(r,rmax,Order);
            j++;
        }
    }   
}

size_t GLCache2D::RAMsize() const
{
    size_t ndoubles=grid.size();
    for (size_t i=0;i<itsDiagGLs_grid_gl.rows();i++)
        for (size_t j=0;j<itsDiagGLs_grid_gl.columns();j++)
        {
            ndoubles+=itsDiagGLs_grid_gl(i,j).RAMsize();
            ndoubles+=itsDiagGLs_gl_grid(j,i).RAMsize();
        }
    return ndoubles;
}


