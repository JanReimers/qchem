// File: IE_Primatives.H get all calculation of primative integrals in one place.
#ifndef _SG_IE_Primatives_H_
#define _SG_IE_Primatives_H_

#include "Imp/BasisSet/AtomIE.H"

namespace SphericalGaussian
{
   
class IE_Primatives 
    : public virtual  Primative_Overlap<double>
    , public virtual  Primative_Kinetic<double>
    , public virtual  Primative_Nuclear<double>
    , public virtual Primative_Repulsion<double>
{
    protected:
    virtual double Overlap(double ea , double eb,size_t l_total) const;
    virtual double Kinetic(double ea , double eb,size_t la, size_t lb) const;
    virtual double Nuclear(double ea , double eb,size_t l_total) const;
    virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const;
};

} //

#endif //_SG_IE_Primatives_H_
