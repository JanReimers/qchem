// File: BasisSet1/Atom/Evaluators/Gaussian/IBS_Evaluator.C
module;
#include <iosfwd>
#include <blaze/Math.h>
export module qchem.BasisSet1.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet1.Atom.Evaluators.Internal.Exponential_IBS_Evaluator;
import qchem.Symmetry.Yl;
import Common.IntPow;

export class Gaussian_IBS_Evaluator : public Exponential_IBS_Evaluator
{
public: 
 
    Gaussian_IBS_Evaluator(const rvec_t& es, int l, const is_t& mls) : Exponential_IBS_Evaluator(es,l,mls) {ns=norms();}
    Gaussian_IBS_Evaluator(const rvec_t& es, int l) : Gaussian_IBS_Evaluator(es,l,{}) {}
    Gaussian_IBS_Evaluator(const rvec_t& es, const Irrep_QNs::sym_t& ir) : Exponential_IBS_Evaluator(es,ir) {ns=norms();}
    Gaussian_IBS_Evaluator(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir) 
    : Exponential_IBS_Evaluator(exponents(N,emin,emax,ir),ir) {ns=norms();}

    Gaussian_IBS_Evaluator Rescale(double scale_factor) const
    {
        return Gaussian_IBS_Evaluator(scale_factor*es,0);
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
    virtual rmat_t XKinetic  (const IBS_Evaluator&) const;

    virtual dERI3  Overlap  (const IBS_Evaluator&) const; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator&) const; //3 center
    
    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;

protected:
    static rvec_t exponents(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir);
    rvec_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const; // RKB needs to override
    static double Grad2(double ea , double eb,size_t la, size_t lb); // RKB needs access
    static double Inv_r2(double ea , double eb,size_t l_total); // RKB needs access
    template <class v> static v gaussian(double r,size_t l,const v& e, const v& n)
    {
        return n*uintpow(r,l)*exp(-e*r*r);
    }
    template <class v> static v grad_gaussian(double r,size_t l,const v& e, const v& n)
    {
        double lr= r==0 ? 0 : l/r;
        return (lr-2*r*e)*gaussian(r,l,e,n);
    }

    

};

export class Gaussian_RKBS_IBS_Evaluator : public Gaussian_IBS_Evaluator
{
public:
    Gaussian_RKBS_IBS_Evaluator(const rvec_t& es, int _kappa, int l, const is_t& mls) : Gaussian_IBS_Evaluator(es,l,mls), kappa(_kappa) {ns=norms();}
    Gaussian_RKBS_IBS_Evaluator(const rvec_t& es, int _kappa, int l) : Gaussian_RKBS_IBS_Evaluator(es,_kappa,l,{}) {}
    Gaussian_RKBS_IBS_Evaluator(size_t N, double emin, double emax, int _kappa, int l): Gaussian_IBS_Evaluator(N,emin,emax,Irrep_QNs::sym_t(new Yl_Sym(0))), kappa(_kappa) {ns=norms();}
    virtual rvec_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;

private:
    rvec_t eval(const rvec3_t&) const;
    int kappa;
};


