// File: SphericalGaussian/IntegralEngine.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/SphericalGaussianCD.H"
#include "Imp/Integrals/AngularIntegrals.H"

namespace SphericalGaussian
{


double Fit_IE::Charge(double ea, size_t l) const
{
    return GaussianIntegral(ea,l);
}

double Orbital_IE::DFTOverlap(double ea, double eb,size_t l) const
{
    return GaussianIntegral(ea+eb,l);
}
double IntegralEngine1::Kinetic(double ea , double eb,size_t l) const
{
    double t=ea+eb;
    size_t l1=l+1;
    return 0.5*(
            (l1*l1 + l*l1) * GaussianIntegral(t,2*l-2)
            -2*l1 * t      * GaussianIntegral(t,2*l  )
            +4*ea*eb       * GaussianIntegral(t,2*l+2)
        );
}

double IntegralEngine1::Nuclear(double ea , double eb,size_t l) const
{
    return GaussianIntegral(ea+eb,2*l-1); //Already has 4*Pi
}
double IntegralEngine1::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    SphericalGaussianCD cd(eab,ec,std::max(la,lc));
    return 4*4*pi*pi*cd.Coulomb_R0(la,lc);
}

double  IntegralEngine1::Overlap(double ea , double eb,size_t l) const
{
    return GaussianIntegral(ea+eb,2*l); //Already has 4*Pi
}

double IntegralEngine1::Integral(qchem::IType type,double ea, double eb,size_t l) const
{
    switch(type)
    {
        case qchem::Overlap1: return GaussianIntegral(ea+eb,2*l); //Already has 4*Pi
        case qchem::Kinetic1:
        {
            double t=ea+eb;
            size_t l1=l+1;
            return 0.5*(
                    (l1*l1 + l*l1) * GaussianIntegral(t,2*l-2)
                    -2*l1 * t      * GaussianIntegral(t,2*l  )
                    +4*ea*eb       * GaussianIntegral(t,2*l+2)
                );
        } 
        case qchem::Nuclear1:  return GaussianIntegral(ea+eb,2*l-1); //Already has 4*Pi
        case qchem::Repulsion1: 
        {
            SphericalGaussianCD cd(ea,eb,std::max(l,l));
            return 4*4*pi*pi*cd.Coulomb_R0(l,l);
        }
        default: assert(false);
    }
    return 0.0;
}


} //namespace
