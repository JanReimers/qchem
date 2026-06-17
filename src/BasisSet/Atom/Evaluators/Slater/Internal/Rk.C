// File: src/BasisSet/Atom/Evaluators/Slater/Internal/Rk.C  Slater integrals for radial Slater functions.
module;
export module qchem.BasisSet.Atom.Evaluators.Slater.Internal.Rk; 
export import qchem.BasisSet.Atom.Evaluators.Internal.Rk;

export namespace Slater
{
//
// Store exponents and derivative tables for a 4 electron charge distribution ab|cd
// formed from a product of slater orbitals.
// eab = ea + eb
// ecd = ec + ed
//
class RkEngine : public virtual Rk
{
public:
    RkEngine(double eab, double ecd, size_t LMax);
    virtual double DirectR0  (size_t la,size_t lc) const; //R_0(la,la,lc,lc);
    virtual double DirectRk  (size_t la,size_t lc,const rvec11_t& Ak) const; //sum{k,A_k*R_k(la,la,lc,lc)};
    virtual double ExchangeRk(size_t la,size_t lb,const rvec11_t& Ak) const; //sum{k,A_k*R_k(la,lb,la,lb)};

    virtual size_t RAMsize() const;
private:
    static double fk(double a, double ab, size_t k,size_t n);
    virtual size_t LMax() const {return itsLMax;}

    double eab, ecd;
    size_t itsLMax;
    rmat_t Iab,Icd; //Derivative tables.
};

} //namespace

