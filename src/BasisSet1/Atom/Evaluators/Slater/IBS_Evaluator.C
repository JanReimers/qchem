// File: BasisSet1/Atom/Evaluators/Slater/IBS_Evaluator.C
module;
#include <iosfwd>
#include <blaze/Math.h>
export module qchem.BasisSet1.Atom.Evaluators.Slater.IBS;
export import qchem.BasisSet1.Atom.Internal.Exponential_IBS_Evaluator;
import Common.IntPow;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Yl;

export class Slater_IBS_Evaluator : public Exponential_IBS_Evaluator
{
public: 
    Slater_IBS_Evaluator(const rvec_t& es, int l, const is_t& mls) : Exponential_IBS_Evaluator(es,l,mls) {ns=norms();}
    Slater_IBS_Evaluator(const rvec_t& es, int l) : Slater_IBS_Evaluator(es,l,{}) {}
    Slater_IBS_Evaluator(const rvec_t& es, const Irrep_QNs::sym_t& ir) : Exponential_IBS_Evaluator(es,ir) {ns=norms();}
    Slater_IBS_Evaluator(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir) 
    : Exponential_IBS_Evaluator(exponents(N,emin,emax,ir),ir) {ns=norms();}
   
    Slater_IBS_Evaluator Rescale(double scale_factor) const
    {
        return Slater_IBS_Evaluator(scale_factor*es,0);
    }
    virtual std::ostream& Write   (std::ostream&) const;
 
    virtual rsmat_t Overlap  () const;
    virtual rsmat_t Grad2    () const;
    virtual rsmat_t Inv_r1   () const;
    virtual rsmat_t Inv_r2   () const;
    virtual rsmat_t Repulsion() const;
    virtual  rvec_t Charge   () const;
    virtual  rvec_t Norm     () const {return ns;}
    virtual rmat_t XRepulsion(const IBS_Evaluator&) const;
    virtual rmat_t XKinetic  (const IBS_Evaluator*) const;

    virtual dERI3  Overlap  (const IBS_Evaluator&) const; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator&) const; //3 center

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;

protected:
    static rvec_t exponents(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir);
    rvec_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const; //Needs an override for RKB.
    static  double Grad2 (double ea , double eb,size_t la, size_t lb); //RKB needs access
    static  double Inv_r2(double ea , double eb,size_t l_total); //RKB needs access

    template <class v> static v slater(double r,size_t l,const v& e, const v& n)
    {
        return n*uintpow(r,l)*exp(-e*r);
    }
    template <class v> static v grad_slater(double r,size_t l,const v& e, const v& n)
    {
        double lr= r==0 ? 0 : l/r;
        return (lr-e)*slater(r,l,e,n);
    }
    


   
};


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