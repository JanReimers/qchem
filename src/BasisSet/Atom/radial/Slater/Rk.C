// File: Slater::RkEngine.H  4 electron Charge distribution of Slater orbitals. 
module;
export module qchem.BasisSet.Atom.Internal.radial.Slater.Rk;
import qchem.BasisSet.Internal.Cache4;
import oml;
export import qchem.Types;

export namespace Slater
{
//
// Store exponents and derivative tables for a 4 electron charge distribution ab|cd
// formed from a product of slater orbitals.
// eab = ea + eb
// ecd = ec + ed
//
class RkEngine : public virtual Cacheable
{
public:
    RkEngine(double eab, double ecd, size_t LMax);
    double Coulomb_R0(size_t la,size_t lc) const; //R_0(la,la,lc,lc);
    Vector<double> Coulomb_Rk(size_t la,size_t lc) const; //R_k(la,la,lc,lc);
    Vector<double> ExchangeRk(size_t la,size_t lb) const; //R_k(la,lb,la,lb);
    //! R_k(la,lb,la,lb) with |Ala-Alb| <= k <= Ala+Alb
    Vector<double> ExchangeRk(size_t Ala,size_t Alb, size_t la,size_t lb) const; 
    
private:
    static double fk(double a, double ab, size_t k,size_t n);
    
    double eab, ecd;
    size_t LMax;
    Matrix<double> Iab,Icd; //Derivative tables.
};

} //namespace

