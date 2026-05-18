// File: BasisSet1/Atom/Evaluators/Gaussian/IBS_Evaluator.C
module;
#include <iosfwd>
#include <blaze/Math.h>
export module qchem.BasisSet1.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet1.Atom.Evaluators.Internal.Exponential_IBS_Evaluator;
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.GaussianIntegrals; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.Rk; 
import qchem.BasisSet1.Atom.Evaluators.Internal.AngularIntegrals;
import qchem.Symmetry.Yl;
import Common.IntPow;

import qchem.BasisSet1.Internal.Cache4;



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
    //using Exponential_IBS_Evaluator::size; does not seem to work for concepts. no member named 'size' in 'Gaussian_BS_Evaluator'
    size_t size() const {return Exponential_IBS_Evaluator::size();}
    virtual std::ostream& Write   (std::ostream&) const;

    double Overlap(size_t i,size_t j) const
    {
        return Gaussian::Integral(es[i]+es[j],2*l)*ns[i]*ns[j]; //Already has 4*Pi and r^2 from dr.
    } 
    double Grad2(size_t i,size_t j) const
    {
        double t=es[i]+es[j];
        size_t l1=l+1;
        return  (l1*l1         * Gaussian::Integral(t,2*l-2)
                -2*l1 * t      * Gaussian::Integral(t,2*l  )
                +4*es[i]*es[j] * Gaussian::Integral(t,2*l+2))*ns[i]*ns[j] ;
    } 
    double Grad2(size_t i,size_t j, const Gaussian_IBS_Evaluator& b) const
    {
        assert(l==b.l);
        double t=es[i]+b.es[j];
        size_t l1=l+1;
        return  (l1*l1           * Gaussian::Integral(t,2*l-2)
                -2*l1 * t        * Gaussian::Integral(t,2*l  )
                +4*es[i]*b.es[j] * Gaussian::Integral(t,2*l+2))*ns[i]*b.ns[j] ;
    } 
    double Inv_r1(size_t i,size_t j) const
    {
        return Gaussian::Integral(es[i]+es[j],2*l-1)*ns[i]*ns[j]; //Already has 4*Pi
    } 
    double Inv_r2(size_t i,size_t j) const
    {
        return Gaussian::Integral(es[i]+es[j],2*l-2)*ns[i]*ns[j]; //Already has 4*Pi
    } 
    double Inv_r2(size_t i,size_t j, const Gaussian_IBS_Evaluator& b) const
    {
        assert(l==b.l);       
        return Gaussian::Integral(es[i]+b.es[j],2*l-2)*ns[i]*b.ns[j]; //Already has 4*Pi
    } 
    double Overlap(size_t i,size_t j, const Gaussian_IBS_Evaluator& c, size_t ic) const
    {
        return Gaussian::Integral(es[i]+es[j]+c.es[ic],2*l+c.l)*ns[i]*ns[j]*c.ns[ic]; //Already has 4*Pi and r^2 from dr.
    } 
    double Repulsion(size_t i,size_t j, const Gaussian_IBS_Evaluator& c, size_t ic) const
    {
        Gaussian::RkEngine cd(es[i]+es[j],c.es[ic],std::max(l,c.l));
        return cd.Coulomb_R0(l,c.l)*FourPi2*ns[i]*ns[j]*c.ns[ic];
    } 
    double Repulsion(size_t i,size_t j) const
    {
        Gaussian::RkEngine cd(es[i],es[j],l);
        return cd.Coulomb_R0(l,l)*FourPi2*ns[i]*ns[j];
    }
    double Repulsion(size_t i,size_t j, const Gaussian_IBS_Evaluator& b) const
    {
        Gaussian::RkEngine cd(es[i],b.es[j],std::max(l,b.l));
        return cd.Coulomb_R0(l,b.l)*FourPi2*ns[i]*b.ns[j];
    }

    double Charge(size_t i) const
    {
        return Gaussian::Integral(es[i],l)*ns[i];
    }
    double Norm(size_t i) const
    {
        return 1.0/sqrt(Gaussian::Integral(2*es[i],2*l));
    }
    double operator()(size_t i, double r) {return gaussian(r,l,es[i],ns[i]);}
    virtual  rvec_t Norm     () const {return ns;}

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;
    virtual std::string RadialType() const;
    virtual Cache41*    MakeCache4() const;
    using Exponential_IBS_Evaluator::maxSpan;
    using rvec11_t=AngularIntegrals::rvec11_t;
    static double direct(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const Gaussian::RkEngine* cd = dynamic_cast<const Gaussian::RkEngine*>(c);
        return cd->Coulomb_Rk(la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const Gaussian::RkEngine* cd = dynamic_cast<const Gaussian::RkEngine*>(c);
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

static_assert(isGeneric_Evaluator<Gaussian_IBS_Evaluator>);
static_assert(is1E_Evaluator     <Gaussian_IBS_Evaluator>);
static_assert(isFit_Evaluator    <Gaussian_IBS_Evaluator>);
static_assert(isDFT_Evaluator    <Gaussian_IBS_Evaluator>);
static_assert(isRKBL_Evaluator   <Gaussian_IBS_Evaluator>);
static_assert(isHF_Evaluator     <Gaussian_IBS_Evaluator>);

export class Gaussian_Cache4 : public  Cache41
{
public:
    // using IBS_Evaluator_t = Gaussian_IBS_Evaluator;
    virtual void Register(Cache4_Client * eval)
    {
        assert(eval);
        Gaussian_IBS_Evaluator* geval=dynamic_cast<Gaussian_IBS_Evaluator*>(eval);
        geval->Register(&grouper);
    }
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const
    {
        return new Gaussian::RkEngine(
            grouper.unique_esv[ia]+grouper.unique_esv[ib],
            grouper.unique_esv[ic]+grouper.unique_esv[id],
            grouper.LMax(ia,ib,ic,id));
    }
private:
    friend class Cache4Tests;
    ExponentGrouper grouper;
};

export class Gaussian_RKBS_IBS_Evaluator : public Gaussian_IBS_Evaluator
{
public:
    Gaussian_RKBS_IBS_Evaluator(const rvec_t& es, int _kappa, int l, const is_t& mls) : Gaussian_IBS_Evaluator(es,l,mls), kappa(_kappa) {ns=norms();}
    Gaussian_RKBS_IBS_Evaluator(const rvec_t& es, int _kappa, int l) : Gaussian_RKBS_IBS_Evaluator(es,_kappa,l,{}) {}
    Gaussian_RKBS_IBS_Evaluator(size_t N, double emin, double emax, int _kappa, int l): Gaussian_IBS_Evaluator(N,emin,emax,Irrep_QNs::sym_t(new Yl_Sym(0))), kappa(_kappa) {ns=norms();}
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

static_assert(isGeneric_Evaluator<Gaussian_RKBS_IBS_Evaluator>);
static_assert(is1E_Evaluator     <Gaussian_RKBS_IBS_Evaluator>);
static_assert(isFit_Evaluator    <Gaussian_RKBS_IBS_Evaluator>);
static_assert(isDFT_Evaluator    <Gaussian_RKBS_IBS_Evaluator>);
static_assert(isRKBL_Evaluator   <Gaussian_RKBS_IBS_Evaluator>);


