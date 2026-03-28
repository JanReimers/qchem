// File: BasisSet/Atom/Gaussian/RKB/IBS_Evaluator.C
module;
export module BasisSet.Atom.Gaussian.RKB.IBS_EValuator;
import BasisSet.Atom.Gaussian.NR.IBS_EValuator;
import qchem.Types;

export class Gaussian_RKBS_IBS : public Gaussian_IBS
{
public:
    Gaussian_RKBS_IBS(const ds_t& es, int _kappa, int l, const is_t& mls) : Gaussian_IBS(es,l,mls), kappa(_kappa) {ns=norms();}
    Gaussian_RKBS_IBS(const ds_t& es, int _kappa, int l) : Gaussian_RKBS_IBS(es,_kappa,l,{}) {}
    virtual ds_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;
private:
    ds_t eval(const rvec3_t&) const;
    int kappa;
};