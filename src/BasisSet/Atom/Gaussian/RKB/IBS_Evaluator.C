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
    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
private:
    ds_t eval(const RVec3&) const;
    int kappa;
};