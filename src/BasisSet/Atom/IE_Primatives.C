// File: src/BasisSet/Atom/IE_Primatives.C  Define primative integrals required for atoms.
export module qchem.BasisSet.Atom.IE_Primatives;
import qchem.Types;

//! Each basis set must implement these integrals.  This version is for exponential, like Slater and Gaussian.
//! ea, eb are the exponents.
export class IE_Primatives
{
public:
    virtual double Overlap  (double ea ,double eb,size_t l_total      ) const=0;
    virtual double Grad2    (double ea ,double eb,size_t la, size_t lb) const=0;
    virtual double Inv_r1   (double ea ,double eb,size_t l_total      ) const=0;
    virtual double Inv_r2   (double ea ,double eb,size_t l_total      ) const=0;
    virtual double Repulsion(double ea ,double ec,size_t la, size_t lc) const=0;
    virtual double Charge   (double ea           ,size_t l            ) const=0;
};
