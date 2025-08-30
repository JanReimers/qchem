// File: BasisSet/Atom/Slater/NR/IBS_Evaluator.C
module;
export module BasisSet.Atom.Slater.RKB.IBS_Evaluator;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;

export class Slater_RKBS_IBS : public Slater_IBS
{
public:
    Slater_RKBS_IBS(const ds_t& es, int _kappa, int l,const is_t& mls) : Slater_IBS(es,l,mls), kappa(_kappa) {ns=norms();}
    Slater_RKBS_IBS(const ds_t& es, int _kappa, int l) : Slater_RKBS_IBS(es,_kappa,l,{}) {}
    ds_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
private:
    ds_t eval(const RVec3&) const;
    int kappa;
};