// File: BasisSet/Atom/Gaussian/NR/IBS_Evaluator.C
module;
#include <valarray>
#include <vector>
#include <iosfwd>
export module BasisSet.Atom.Gaussian.NR.IBS_EValuator;
import qchem.BasisSet.Atom.Internal.Exponential_IBS_Evaluator;
import Common.IntPow;

export class Gaussian_IBS : public Exponential_IBS_Evaluator
{
public: 
 
    Gaussian_IBS(const   ds_t& es, int l, const is_t& mls) : Exponential_IBS_Evaluator(es,l,mls) {ns=norms();}
    Gaussian_IBS(const   ds_t& es, int l) : Gaussian_IBS(es,l,{}) {}

    virtual std::ostream& Write   (std::ostream&) const;

    virtual omls_t Overlap  () const;
    virtual omls_t Grad2    () const;
    virtual omls_t Inv_r1   () const;
    virtual omls_t Inv_r2   () const;
    virtual omls_t Repulsion() const;
    virtual omlv_t Charge   () const;
    virtual ds_t   Norm     () const {return ns;}
    virtual omlm_t XRepulsion(const Fit_IBS&) const;
    virtual omlm_t XKinetic  (const Orbital_RKBS_IBS<double>*) const;

    virtual dERI3  Overlap  (const Fit_IBS&) const; //3 center
    virtual dERI3  Repulsion(const Fit_IBS&) const; //3 center
    
    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

protected:
    ds_t norms() const; //assumes es,l are already initialized
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

    template <class T> static Vector<T> convert(const std::valarray<T>& v) 
    {
        Vector<T> ret(v.size());
        size_t i=0;
        for (auto vi:v) ret(++i)=vi;
        return ret;
    }

};


