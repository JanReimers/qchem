// File: BasisSet/Atom/Slater/NR/IBS_Evaluator.C
module;
#include <string>
export module BasisSet.Atom.Slater.RKB.IBS_Evaluator;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Yl;

export class Slater_RKBS_IBS_Evaluator : public Slater_IBS_Evaluator
{
public:
    Slater_RKBS_IBS_Evaluator(const rvec_t& es, int _kappa, int l,const is_t& mls) : Slater_IBS_Evaluator(es,l,mls), kappa(_kappa) {ns=norms();}
    Slater_RKBS_IBS_Evaluator(const rvec_t& es, int _kappa, int l) : Slater_RKBS_IBS_Evaluator(es,_kappa,l,{}) {}
    Slater_RKBS_IBS_Evaluator(size_t N, double emin, double emax, int _kappa, int l): Slater_IBS_Evaluator(N,emin,emax,Irrep_QNs::sym_t(new Yl_Sym(0))), kappa(_kappa) {ns=norms();}
    rvec_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;
private:
    rvec_t eval(const rvec3_t&) const;
    int kappa;
};