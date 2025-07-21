// File: Slater/IE_Primatives.H get all calculation of primative integrals in one place.
#ifndef _Slater_IE_Primatives_H_
#define _Slater_IE_Primatives_H_

import qchem.BasisSet.Atom.IE;

namespace Slater
{
    class IE_Primatives 
    : public virtual Primative_Overlap  <double>
    , public virtual Primative_Grad2    <double>
    , public virtual Primative_Inv_r1   <double>
    , public virtual Primative_Inv_r2   <double>
    , public virtual Primative_Repulsion<double>
    , public virtual Primative_Charge   <double>
{
    protected:
    virtual double Overlap  (double ea, double eb,size_t l_total) const;
    virtual double Grad2    (double ea, double eb,size_t la, size_t lb) const;
    virtual double Inv_r1   (double ea, double eb,size_t l_total) const;
    virtual double Inv_r2   (double ea, double eb,size_t l_total) const;
    virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const;
    virtual double Charge   (double ea,           size_t l) const;
};

}

#endif // _Slater_IE_Primatives_H_
   