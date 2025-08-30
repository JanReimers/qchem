// File: BasisSet/Atom/Rk.C Interface for all Rk (Slater Integral) engines.
module;

export module qchem.BasisSet.Atom.Rk;
export import qchem.BasisSet.Internal.Cache4;
export import qchem.Types;
//
//  These are often called Slater integrals. They represent the radial part of the 
//  2 electron repulsion integrals (ERIs) encountered in atomic Hartree-Fock (HF) 
//  calculations.  For large l they can become rather complicated. But recursion
//  relations allows us to calculate and store tables of these integrals in an 
//  efficient manner.
//  Using this interface is part of the architecture for enabling a generic
//  quadruple nested loop required for evaluating direct and exchange integrals
//  over gaussian/slater/b-spline irrep basis sets.
//
export class Rk : public virtual Cacheable
{
public:
    virtual ~Rk() {};
    virtual double Coulomb_R0(size_t la,size_t lc) const=0; //R_0(la,la,lc,lc);
    virtual RVec   Coulomb_Rk(size_t la,size_t lc) const=0; //R_k(la,la,lc,lc);
    virtual RVec   ExchangeRk(size_t la,size_t lb) const=0; //R_k(la,lb,la,lb);
};