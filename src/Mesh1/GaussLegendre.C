// File: GaussLegendre.C  Shared 1D Gauss-Legendre quadrature.
//
// ONE implementation, to be used by the radial meshes, the angular Gauss-Legendre scheme, the PW
// real-space grid, and (on migration) the Atom BSpline integrator -- each of which historically
// rolled its own.  Ported from gauleg.f (Numerical Recipes), the canonical source.
//
// NB: the old src/Mesh angular copy wrote x[n+1-i+1] (== x[n+1] for i=1), an out-of-bounds index
// on a size-n vector.  The Fortran uses x(n+1-i), i.e. 0-based x[n-i].  We follow the Fortran.
module;
#include <cmath>
export module qchem.Mesh1.GaussLegendre;
export import qchem.Types;
import qchem.Math;

//! \brief n-point Gauss-Legendre nodes \c x and weights \c w on [a,b].  Exact for polynomials
//! up to degree 2n-1.
export struct GaussLegendre
{
    rvec_t x, w;
    GaussLegendre(int n, double a, double b);
};

GaussLegendre::GaussLegendre(int n, double a, double b)
    : x(n), w(n)
{
    const double eps=3.0e-14;
    int    m =(n+1)/2;                 // roots are symmetric -> compute half
    double xm=0.5*(b+a), xl=0.5*(b-a);
    for (int i=1; i<=m; i++)
    {
        double z=cos(Pi*(i-0.25)/(n+0.5)), z1, pp;
        do
        {
            double p1=1.0, p2=0.0;
            for (int j=1; j<=n; j++)   // Legendre recurrence -> P_n(z)
            {
                double p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.0);  // P_n'(z)
            z1=z;
            z=z1-p1/pp;                // Newton step
        }
        while (std::fabs(z-z1) > eps);
        x[i-1]  = xm-xl*z;             // symmetric pair, 0-based (cf. Fortran x(i), x(n+1-i))
        x[n-i]  = xm+xl*z;
        w[i-1]  = 2.0*xl/((1.0-z*z)*pp*pp);
        w[n-i]  = w[i-1];
    }
}
