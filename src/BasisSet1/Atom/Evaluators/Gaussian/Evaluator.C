// File: BasisSet1/Atom/Evaluators/Gaussian/Evaluator.C
module;
#include <iosfwd>
#include <blaze/Math.h>
export module qchem.BasisSet.Atom.Evaluators.Gaussian.IBS;
import qchem.BasisSet.Atom.Evaluators.Internal.ExponentialEvaluator;
import qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
import qchem.BasisSet.Atom.Evaluators.Internal.RKBL_Angular;
import qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.GaussianIntegrals;
import qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.Concepts;
import qchem.Symmetry.Spherical;
import qchem.IntPow;
import qchem.BasisSet.Internal.Cache4;


export namespace BasisSet::Atom::Evaluators::Gaussian
{

// Gaussian::Radial holds all Gaussian-specific 1e integrals and Rk machinery,
// shared between NR (Evaluator) and RKB large component (RKBL_Evaluator).
class Radial : public ExponentialEvaluator
{
public:
    Radial(const rvec_t& _es, int _l) 
    : ExponentialEvaluator(_es,_l)
    , l(_l) 
    {ns=norms();}
    Radial(const rvec_t& _es, const sym_t& ir, size_t ltrim=0)
        : ExponentialEvaluator(_es,ir,ltrim) 
        , l(Symmetry::Getl(ir))
        {
            ns=norms();
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
    double Grad2(size_t i,size_t j, const Radial& s) const
    {
        assert(l==s.l);
        double t=es[i]+s.es[j];
        size_t l1=l+1;
        return  (l1*l1           * ::Gaussian::Integral(t,2*l-2)
                -2*l1 * t        * ::Gaussian::Integral(t,2*l  )
                +4*es[i]*s.es[j] * ::Gaussian::Integral(t,2*l+2))*ns[i]*s.ns[j] ;
    }
    double Inv_r1(size_t i,size_t j) const
    {
        return ::Gaussian::Integral(es[i]+es[j],2*l-1)*ns[i]*ns[j]; //Already has 4*Pi
    }
    double Inv_r2(size_t i,size_t j) const
    {
        return ::Gaussian::Integral(es[i]+es[j],2*l-2)*ns[i]*ns[j]; //Already has 4*Pi
    }
    double Inv_r2(size_t i,size_t j, const Radial& b) const
    {
        assert(l==b.l);
        return ::Gaussian::Integral(es[i]+b.es[j],2*l-2)*ns[i]*b.ns[j]; //Already has 4*Pi
    }
    double Overlap(size_t i,size_t j, const Radial& c, size_t ic) const
    {
        return ::Gaussian::Integral(es[i]+es[j]+c.es[ic],2*l+c.l)*ns[i]*ns[j]*c.ns[ic]; //Already has 4*Pi and r^2 from dr.
    }
    double Repulsion(size_t i,size_t j, const Radial& c, size_t ic) const
    {
        ::Gaussian::RkEngine cd(es[i]+es[j],c.es[ic],std::max(l,c.l));
        return cd.Coulomb_R0(l,c.l)*FourPi2*ns[i]*ns[j]*c.ns[ic];
    }
    double Repulsion(size_t i,size_t j) const
    {
        ::Gaussian::RkEngine cd(es[i],es[j],l);
        return cd.Coulomb_R0(l,l)*FourPi2*ns[i]*ns[j];
    }
    double Repulsion(size_t i,size_t j, const Radial& b) const
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
    static rvec_t exponents(size_t N, double emin, double emax, const sym_t& ir);
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

    int l;
};

// NR HF evaluator: Gaussian radial + NR angular.
class Evaluator : public Radial, public NR_Angular
{
public:
    // Used only for the rescale operation to make Fit basis sets.
    Evaluator(const rvec_t& es, int l, const ivec_t& mls={})
        : Radial(es,l), NR_Angular(l,mls) {}

        Evaluator(const rvec_t& es, const sym_t& ir, size_t ltrim=0)
        : Radial(es,ir,ltrim)
        , NR_Angular(ir) {}
    Evaluator(size_t N, double emin, double emax, const sym_t& ir)
        : Radial(Radial::exponents(N,emin,emax,ir),ir)
        , NR_Angular(ir) {}

    Evaluator Rescale(double scale_factor) const { return Evaluator(scale_factor*es,Getl()); }
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

class RKBL_Evaluator : public RKB_Angular, public Radial
{
public:
    // Used only for the rescale operation to make Fit basis sets.
    // RKBL_Evaluator(const rvec_t& es, int l, const ivec_t& mls={})
    //     : Radial(es,l), RKB_Angular(l,mls) {}

    // RKBL_Evaluator(const rvec_t& es, const sym_t& ir, size_t ltrim=0)
    //     : Radial(es,ir,ltrim)
    //     , RKB_Angular(ir) {}
    RKBL_Evaluator(size_t N, double emin, double emax, const sym_t& ir)
        : RKB_Angular(ir)
        , Radial(Radial::exponents(N,emin,emax,ir),ir)
        {}

    // Evaluator Rescale(double scale_factor) const { return Evaluator(scale_factor*es,Getl()); }
};

// RKB small-component evaluator: shares Gaussian radial with NR Evaluator.
class RKBS_Evaluator : public RKBL_Evaluator
{
public:
    // RKBS_Evaluator(const rvec_t& es, int _κ, int l)
    //     : RKB_Angular(ir)
    //     , Radial(es,ir), κ(Getκ()) 
    //     {ns=norms();}
    RKBS_Evaluator(size_t N, double emin, double emax, const sym_t& ir)
        : RKBL_Evaluator(N,emin,emax,ir)
        {ns=norms();}
    // RKBS_Evaluator(const rvec_t& es, const sym_t& ir)
    //     : RKB_Angular(ir)
    //     , Radial(es,ir)
    //     {ns=norms();}
    // int Getκ() const { return κ; }
    // RKBS_Evaluator(size_t N, double emin, double emax, int κ, int l);
    virtual rvec_t norms() const;
    double Inv_r1(size_t i,size_t j) const
    {
        double p=es[i]+es[j];
        return 8*Pi*es[i]*es[j]/(p*p)*ns[i]*ns[j];
    }

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;

private:
    using Radial::l;
    // int    κ;
    rvec_t eval(const rvec3_t&) const;
};

static_assert(isGeneric_Evaluator<RKBS_Evaluator>);
static_assert(is1E_Evaluator     <RKBS_Evaluator>);
static_assert(isFit_Evaluator    <RKBS_Evaluator>);
static_assert(isDFT_Evaluator    <RKBS_Evaluator>);
static_assert(isRKBL_Evaluator   <RKBS_Evaluator>);


} //namespace