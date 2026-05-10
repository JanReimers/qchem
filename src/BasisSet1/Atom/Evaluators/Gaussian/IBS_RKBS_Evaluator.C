// File: BasisSet/Atom/Gaussian/RKB/IBS_Evaluator.C
module;
#include <string>
export module BasisSet.Atom.Gaussian.RKB.IBS_EValuator;
import BasisSet.Atom.Gaussian.NR.IBS_EValuator;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Yl;
import qchem.Types;

export class Gaussian_RKBS_IBS : public Gaussian_IBS
{
public:
    Gaussian_RKBS_IBS(const rvec_t& es, int _kappa, int l, const is_t& mls) : Gaussian_IBS(es,l,mls), kappa(_kappa) {ns=norms();}
    Gaussian_RKBS_IBS(const rvec_t& es, int _kappa, int l) : Gaussian_RKBS_IBS(es,_kappa,l,{}) {}
    Gaussian_RKBS_IBS(size_t N, double emin, double emax, int _kappa, int l): Gaussian_IBS(N,emin,emax,Irrep_QNs::sym_t(new Yl_Sym(0))), kappa(_kappa) {ns=norms();}
    virtual rvec_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;

private:
    rvec_t eval(const rvec3_t&) const;
    int kappa;
};