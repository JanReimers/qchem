// File: BasisSet1/Atom/Evaluators/Gaussian/Evaluator.C
module;
#include <iosfwd>
#include <blaze/Math.h>
export module qchem.BasisSet.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Atom.Evaluators.Internal.ExponentialEvaluator;
import qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.GaussianIntegrals; 
import qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.Rk; 
import qchem.BasisSet.Atom.Evaluators.Concepts;
import qchem.Symmetry.Yl;
import qchem.IntPow;

import qchem.BasisSet.Internal.Cache4;


export namespace BasisSet::Atom::Evaluators::Gaussian
{
class Evaluator : public ExponentialEvaluator
{
public: 
 
    Evaluator(const rvec_t& es, int l, const ivec_t& mls) : ExponentialEvaluator(es,l,mls) {ns=norms();}
    Evaluator(const rvec_t& es, int l) : Evaluator(es,l,{}) {}
    Evaluator(const rvec_t& es, const Irrep_QNs::sym_t& ir, size_t ltrim=0) : ExponentialEvaluator(es,ir,ltrim) {ns=norms();}
    Evaluator(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir) 
    : ExponentialEvaluator(exponents(N,emin,emax,ir),ir) {ns=norms();}

    Evaluator Rescale(double scale_factor) const
    {
        return Evaluator(scale_factor*es,0);
    }
    virtual std::ostream& Write   (std::ostream&) const;

    double Overlap(size_t i,size_t j) const
    {
        return ::Gaussian::Integral(es[i]+es[j],2*l)*ns[i]*ns[j]; //Already has 4*Pi and r^2 from dr.
    } 
    double Grad2(size_t i,size_t j) const
    {
        double t=es[i]+es[j];
        size_t l1=l+1;
        return  (l1*l1         * ::Gaussian::Integral(t,2*l-2)
                -2*l1 * t      * ::Gaussian::Integral(t,2*l  )
                +4*es[i]*es[j] * ::Gaussian::Integral(t,2*l+2))*ns[i]*ns[j] ;
    } 
    double Grad2(size_t i,size_t j, const Evaluator& b) const
    {
        assert(l==b.l);
        double t=es[i]+b.es[j];
        size_t l1=l+1;
        return  (l1*l1           * ::Gaussian::Integral(t,2*l-2)
                -2*l1 * t        * ::Gaussian::Integral(t,2*l  )
                +4*es[i]*b.es[j] * ::Gaussian::Integral(t,2*l+2))*ns[i]*b.ns[j] ;
    } 
    double Inv_r1(size_t i,size_t j) const
    {
        return ::Gaussian::Integral(es[i]+es[j],2*l-1)*ns[i]*ns[j]; //Already has 4*Pi
    } 
    double Inv_r2(size_t i,size_t j) const
    {
        return ::Gaussian::Integral(es[i]+es[j],2*l-2)*ns[i]*ns[j]; //Already has 4*Pi
    } 
    double Inv_r2(size_t i,size_t j, const Evaluator& b) const
    {
        assert(l==b.l);       
        return ::Gaussian::Integral(es[i]+b.es[j],2*l-2)*ns[i]*b.ns[j]; //Already has 4*Pi
    } 
    double Overlap(size_t i,size_t j, const Evaluator& c, size_t ic) const
    {
        return ::Gaussian::Integral(es[i]+es[j]+c.es[ic],2*l+c.l)*ns[i]*ns[j]*c.ns[ic]; //Already has 4*Pi and r^2 from dr.
    } 
    double Repulsion(size_t i,size_t j, const Evaluator& c, size_t ic) const
    {
        ::Gaussian::RkEngine cd(es[i]+es[j],c.es[ic],std::max(l,c.l));
        return cd.Coulomb_R0(l,c.l)*FourPi2*ns[i]*ns[j]*c.ns[ic];
    } 
    double Repulsion(size_t i,size_t j) const
    {
        ::Gaussian::RkEngine cd(es[i],es[j],l);
        return cd.Coulomb_R0(l,l)*FourPi2*ns[i]*ns[j];
    }
    double Repulsion(size_t i,size_t j, const Evaluator& b) const
    {
        ::Gaussian::RkEngine cd(es[i],b.es[j],std::max(l,b.l));
        return cd.Coulomb_R0(l,b.l)*FourPi2*ns[i]*b.ns[j];
    }

    double Charge(size_t i) const
    {
        return ::Gaussian::Integral(es[i],l)*ns[i];
    }
    double Norm(size_t i) const
    {
        return 1.0/sqrt(::Gaussian::Integral(2*es[i],2*l));
    }
    // double operator()(size_t i, double r) {return gaussian(r,l,es[i],ns[i]);}
    virtual  rvec_t Norm     () const {return ns;}

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;
    virtual std::string RadialType() const;
    virtual Cache4*    MakeCache4() const;
    using rvec11_t=rvec11_t;
    static double direct(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::Gaussian::RkEngine* cd = dynamic_cast<const ::Gaussian::RkEngine*>(c);
        return cd->Coulomb_Rk(la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::Gaussian::RkEngine* cd = dynamic_cast<const ::Gaussian::RkEngine*>(c);
        return cd->ExchangeRk(la,lc,Ak); // contract over k Rk*Ak, exchange version is more complicated
    }

protected:
    static rvec_t exponents(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir);
    rvec_t norms() const; //assumes es,l are already initialized
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

static_assert(isGeneric_Evaluator<Evaluator>);
static_assert(is1E_Evaluator     <Evaluator>);
static_assert(isFit_Evaluator    <Evaluator>);
static_assert(isDFT_Evaluator    <Evaluator>);
static_assert(isRKBL_Evaluator   <Evaluator>);
static_assert(isHF_Evaluator     <Evaluator>);

class Gaussian_Cache4 : public  Cache4
{
public:
    // using Evaluator_t = Evaluator;
    virtual void Register(Cache4_Client * eval)
    {
        assert(eval);
        Evaluator* geval=dynamic_cast<Evaluator*>(eval);
        geval->Register(&grouper); //Look for new exponents and get indies for exponents in geval.
        //
        //  At this point we need sweep through all Cacheable* (Rks) in Cache4::cache_t
        //  and check if geval is supported (geval.l <= Rk.LMax).
        //  All unsupport Rks will be removed.  These will then automatically be recreated next time
        //  loop_4 is called.
        //
        ::Cache4::Register(eval);
    }
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const
    {
        return new ::Gaussian::RkEngine(
            grouper.unique_esv[ia]+grouper.unique_esv[ib],
            grouper.unique_esv[ic]+grouper.unique_esv[id],
            grouper.LMax(ia,ib,ic,id));
    }
private:
    friend class Cache4Tests;
    ExponentGrouper grouper;
};

class RKBS_Evaluator : public Evaluator
{
public:
    RKBS_Evaluator(const rvec_t& es, int _kappa, int l, const ivec_t& mls) : Evaluator(es,l,mls), kappa(_kappa) {ns=norms();}
    RKBS_Evaluator(const rvec_t& es, int _kappa, int l) : RKBS_Evaluator(es,_kappa,l,{}) {}
    RKBS_Evaluator(size_t N, double emin, double emax, int _kappa, int l): Evaluator(N,emin,emax,Irrep_QNs::sym_t(new Yl_Sym(0))), kappa(_kappa) {ns=norms();}
    virtual rvec_t norms() const; //assumes es,l are already initialized
    double Inv_r1(size_t i,size_t j) const
    {
        return 4*es[i]*es[j]*::Gaussian::Integral(es[i]+es[j],2*l+1)*ns[i]*ns[j]; //Already has 4*Pi
    }    
    

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;

private:
    rvec_t eval(const rvec3_t&) const;
    int kappa;
};

static_assert(isGeneric_Evaluator<RKBS_Evaluator>);
static_assert(is1E_Evaluator     <RKBS_Evaluator>);
static_assert(isFit_Evaluator    <RKBS_Evaluator>); 
static_assert(isDFT_Evaluator    <RKBS_Evaluator>); 
static_assert(isRKBL_Evaluator   <RKBS_Evaluator>);


} //namespace