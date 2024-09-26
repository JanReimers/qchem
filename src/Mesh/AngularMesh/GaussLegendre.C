#include "oml/vector.h"
#include "Misc/DFTDefines.H"
#include <cmath>


void GaussLegendre(double x1, double x2, Vector<double>& x, Vector<double>& w, int n)
{
    const double eps=3.0E-14;
    double pp=0, z=0, z1=0;
    int m=(n+1)/2;
    double xm=0.5*(x2+x1), xl=0.5*(x2-x1);
    for (int i=1; i<=m; i++)
    {
        z=cos(Pi*(i-0.25)/(n+0.5));
        do
        {
            double p1=1.0,p2=0.0,p3=0.0;
            for (int j=1; j<=n; j++)
            {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        }
        while (fabs(z-z1) > eps);
        x(i    )=xm-xl*z;
        x(n+1-i)=xm+xl*z;
        w(i)=2.0*xl/((1.0-z*z)*pp*pp);
        w(n+1-i)=w(i);
    }
}
