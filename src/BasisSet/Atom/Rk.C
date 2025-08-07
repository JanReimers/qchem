// File: BasisSet/Atom/Rk.C Interface for all Rk engines.
module;

export module qchem.BasisSet.Atom.internal.Rk;
export import qchem.Types;
//
//  These are often called Slater integrals. They represent the radial part of the 
//  2 electron repulsion integrals (ERIs) encountered in atomic Hartree-Fock (HF) 
//  calculations.  For large l they can become very complicated. But recursion
//  relations allows to calculate and store tables of the integrals in and 
//  efficient manner.
//  Using this interface is part of the architecture for anabling a generic
//  quadruple nested loop required for evaluating direct and exchange integrals
//  over gaussian/slater/b-spline irrep basis sets.
//
export class Rk
{
public:
    virtual ~Rk() {};
    virtual        double  Coulomb_R0(size_t la,size_t lc) const=0; //R_0(la,la,lc,lc);
    virtual Vector<double> Coulomb_Rk(size_t la,size_t lc) const=0; //R_k(la,la,lc,lc);
    virtual Vector<double> ExchangeRk(size_t la,size_t lb) const=0; //R_k(la,lb,la,lb);
    //! R_k(la,lb,la,lb) with |Ala-Alb| <= k <= Ala+Alb
    //virtual Vector<double> ExchangeRk(size_t Ala,size_t Alb, size_t la,size_t lb) const=0;  is this used???
};