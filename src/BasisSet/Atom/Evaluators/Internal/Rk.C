// File: BasisSet/Atom/Evaluators/Internal/Rk.C Interface for all Rk Slater Integrals.
module;

export module qchem.BasisSet.Atom.Evaluators.Internal.Rk;
export import qchem.BasisSet.Internal.Cache4;
export import qchem.Types;

namespace qchem {
//
//  These are often called Slater integrals. They represent the radial part of the 
//  2 electron repulsion integrals (ERIs) encountered in atomic Hartree-Fock (HF) 
//  calculations.  For large l they can become rather complicated. But recursion
//  relations allow us to calculate and store tables of these integrals in an 
//  efficient manner.
//  Using this interface is part of the architecture for enabling a generic
//  quadruple nested loop required for evaluating direct and exchange integrals
//  over gaussian/slater/b-spline irrep basis sets.
//
export class Rk : public virtual Cacheable4
{
public:
    typedef rvec11_t rvec11_t;
    virtual ~Rk() {};
    virtual bool   isSupported(const Cache4_Client*) const;
    virtual double DirectR0  (size_t la,size_t lc) const=0; //R_0(la,la,lc,lc);
    virtual double DirectRk  (size_t la,size_t lc,const rvec11_t& Ak) const=0; //sum{k,A_k*R_k(la,la,lc,lc)};
    virtual double ExchangeRk(size_t la,size_t lb,const rvec11_t& Ak) const=0; //sum{k,A_k*R_k(la,lb,la,lb)};
private:
    virtual size_t LMax() const=0;
};
} // namespace qchem