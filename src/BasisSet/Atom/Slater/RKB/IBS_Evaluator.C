// File: BasisSet/Atom/Slater/NR/IBS_Evaluator.C
module;
export module BasisSet.Atom.Slater.RKB.IBS_Evaluator;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;

export class Slater_RKBS_IBS : public Slater_IBS
{
public:
    Slater_RKBS_IBS(const rvec_t& es, int _kappa, int l,const is_t& mls) : Slater_IBS(es,l,mls), kappa(_kappa) {ns=norms();}
    Slater_RKBS_IBS(const rvec_t& es, int _kappa, int l) : Slater_RKBS_IBS(es,_kappa,l,{}) {}
    rvec_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;
private:
    rvec_t eval(const rvec3_t&) const;
    int kappa;
};