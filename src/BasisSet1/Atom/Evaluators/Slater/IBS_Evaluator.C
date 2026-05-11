// File: BasisSet/Atom/Slater/NR/IBS_Evaluator.C
module;
#include <iosfwd>
#include <blaze/Math.h>
export module BasisSet.Atom.Slater.NR.IBS_Evaluator;
export import qchem.BasisSet1.Atom.Internal.Exponential_IBS_Evaluator;
import Common.IntPow;

export class Slater_IBS : public Exponential_IBS_Evaluator
{
public: 
    Slater_IBS(const rvec_t& es, int l, const is_t& mls) : Exponential_IBS_Evaluator(es,l,mls) {ns=norms();}
    Slater_IBS(const rvec_t& es, int l) : Slater_IBS(es,l,{}) {}
    Slater_IBS(const rvec_t& es, const Irrep_QNs::sym_t& ir) : Exponential_IBS_Evaluator(es,ir) {ns=norms();}
    Slater_IBS(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir) 
    : Exponential_IBS_Evaluator(exponents(N,emin,emax,ir),ir) {ns=norms();}
   
    Slater_IBS Rescale(double scale_factor) const
    {
        return Slater_IBS(scale_factor*es,0);
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

