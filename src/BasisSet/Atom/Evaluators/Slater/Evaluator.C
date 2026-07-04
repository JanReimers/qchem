// File: BasisSet/Atom/Evaluators/Slater/Evaluator.C
module;
#include <iosfwd>
#include <cassert>
#include "forward.H"
export module qchem.BasisSet.Atom.Evaluators.Slater.IBS;
export import qchem.BasisSet.Atom.Evaluators.Internal.ExponentialEvaluator;
import qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
import qchem.BasisSet.Atom.Evaluators.Internal.RKBL_Angular;
import qchem.BasisSet.Atom.Evaluators.Slater.Internal.Integrals;
import qchem.BasisSet.Atom.Evaluators.Slater.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.Internal.Grouper;
import qchem.BasisSet.Atom.Evaluators;
import qchem.BasisSet.Internal.Cache4;
import qchem.Symmetry.Atom.Spherical;
import qchem.IntPow;
import qchem.Blaze;
import qchem.Math;

export namespace qchem::BasisSet::Atom::Evaluators::Slater
{

// Slater::Radial holds all Slater-specific 1e integrals and Rk machinery,
// shared between NR (Evaluator) and RKB large component (RKBS_Evaluator).
class Radial : public ExponentialEvaluator
{
public:
    Radial(const rvec_t& _es, const sym_t& ir, size_t ltrim=0)
        : ExponentialEvaluator(_es,ir,ltrim) 
        , l(Symmetry::Getl(ir))
        {ns=norms();}
    Radial(size_t N, double emin, double emax, const sym_t& ir, size_t ltrim=0)
        : ExponentialEvaluator(exponents(N,emin,emax,ir),ir,ltrim) 
        , l(Symmetry::Getl(ir))
        {ns=norms();}

    virtual std::ostream& Write   (std::ostream&) const;
    double Overlap(size_t i,size_t j) const
    {
        return ::qchem::Slater::Integral(es[i]+es[j],2*l)*ns[i]*ns[j]; //Already has 4*Pi and r^2 from dr.
    } 
    double Inv_r1(size_t i,size_t j) const
    {
        return ::qchem::Slater::Integral(es[i]+es[j],2*l-1)*ns[i]*ns[j]; //Already has 4*Pi
    }
    double Overlap(size_t i,size_t j, const Radial& c, size_t ic) const
    {
        return ::qchem::Slater::Integral(es[i]+es[j]+c.es[ic],2*l+c.l)*ns[i]*ns[j]*c.ns[ic]; //Already has 4*Pi and r^2 from dr.
    }

        double Grad2(size_t i,size_t j) const
    {
        double ab=es[i]+es[j];
        double Term1=(l+1)*(l+1)*::qchem::Slater::Integral(ab,2*l-2); //SlaterIntegral already has 4*Pi
        double Term2=-(l+1)*ab  *::qchem::Slater::Integral(ab,2*l-1);
        double Term3=es[i]*es[j]*::qchem::Slater::Integral(ab,2*l);
        return (Term1+Term2+Term3)*ns[i]*ns[j];
    } 
    double Inv_r2(size_t i,size_t j) const
    {
        return ::qchem::Slater::Integral(es[i]+es[j],2*l-2)*ns[i]*ns[j]; //Already has 4*Pi
    }

    double Repulsion(size_t i,size_t j, const Radial& c, size_t ic) const
    {
        ::qchem::Slater::RkEngine cd(es[i]+es[j],c.es[ic],max(l,c.l));
        return cd.DirectR0  (l,c.l)*FourPi2*ns[i]*ns[j]*c.ns[ic];
    }
    double Repulsion(size_t i,size_t j) const
    {
        ::qchem::Slater::RkEngine cd(es[i],es[j],l);
        return cd.DirectR0  (l,l)*FourPi2*ns[i]*ns[j];
    }
    double Repulsion(size_t i,size_t j, const Radial& b) const
    {
        ::qchem::Slater::RkEngine cd(es[i],b.es[j],max(l,b.l));
        return cd.DirectR0  (l,b.l)*FourPi2*ns[i]*b.ns[j];
    }
    double Charge(size_t i) const
    {
        return ::qchem::Slater::Integral(es[i],l)*ns[i];
    }
    double Norm(size_t i) const
    {
        return 1.0/sqrt(::qchem::Slater::Integral(2*es[i],2*l));
    }

    virtual  rvec_t Norm     () const {return ns;}

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;
    virtual std::string RadialType() const;
    virtual Cache4*    MakeCache4() const;
    using rvec11_t=rvec11_t;
    static double direct(const Cacheable4* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::qchem::Slater::RkEngine* cd = dynamic_cast<const ::qchem::Slater::RkEngine*>(c);
        return cd->DirectRk  (la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable4* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::qchem::Slater::RkEngine* cd = dynamic_cast<const ::qchem::Slater::RkEngine*>(c);
        return cd->ExchangeRk(la,lc,Ak); // contract over k Rk*Ak, exchange version is more complicated
    }

    int l;

protected:
    rvec_t norms() const; //assumes es,l are already initialized

    template <class v> static v slater(double r,size_t l,const v& e, const v& n)
    {
        return n*uintpow(r,l)*blazem::exp(-e*r);
    }
    template <class v> static v grad_slater(double r,size_t l,const v& e, const v& n)
    {
        double lr= r==0 ? 0 : l/r;
        return (lr-e)*slater(r,l,e,n);
    }

private:
    static rvec_t exponents(size_t N, double emin, double emax, const sym_t& ir);
};

// NR HF evaluator: Slater radial + NR angular.
class NR_Evaluator : public Radial, public NR_Angular
{
public:
    NR_Evaluator(const rvec_t& es, const sym_t& ir, size_t ltrim=0)
        : Radial(es,ir,ltrim)
        , NR_Angular(ir) {}
    NR_Evaluator(size_t N, double emin, double emax, const sym_t& ir, size_t ltrim=0)
        : Radial(N,emin,emax,ir,ltrim)
        , NR_Angular(ir) {}

    NR_Evaluator Rescale(double scale_factor, sym_t s) const { return NR_Evaluator(scale_factor*es,s); }
    virtual int Getl() const override {return NR_Angular::Getl();}
    using Radial::l;
};

static_assert(isOpr_Evaluator<NR_Evaluator>);
static_assert( is1E_Evaluator<NR_Evaluator>);
static_assert(isFit_Evaluator<NR_Evaluator>);
static_assert(isDFT_Evaluator<NR_Evaluator>);
static_assert( isHF_Evaluator<NR_Evaluator>);



// RKB small-component evaluator: shares Slater radial with NR Evaluator.

class RKBL_Evaluator : public RKB_Angular, public Radial
{
public:
    RKBL_Evaluator(size_t N, double emin, double emax, const sym_t& ir)
        : RKB_Angular(ir)
        , Radial(N,emin,emax,ir)
        {}
    RKBL_Evaluator(const rvec_t& es, const sym_t& ir, size_t ltrim=0)
        : RKB_Angular(ir)
        , Radial(es,ir,ltrim)
        {}
    
    virtual int Getl() const override {return RKB_Angular::Getl();}
};

static_assert(isOpr_Evaluator<RKBL_Evaluator>);
static_assert( is1E_Evaluator<RKBL_Evaluator>);
static_assert( isHF_Evaluator<RKBL_Evaluator>);

class RKBS_Evaluator : public RKBL_Evaluator
{
public:
    RKBS_Evaluator(size_t N, double emin, double emax, const sym_t& ir)
        : RKBL_Evaluator(N,emin,emax,ir)
        {ns=norms();}
    RKBS_Evaluator(const rvec_t& es, const sym_t& ir, size_t ltrim=0)
        : RKBL_Evaluator(es,ir,ltrim)
        {ns=norms();}
    rvec_t norms() const;

    double Inv_r1(size_t i,size_t j) const
    {
        // Small-component nuclear attraction <Q|1/r|Q> with Q=((l+1+κ)/r - e)r^l e^-er.
        // The κ-dependent terms (the spin-orbit piece) vanish for j=l+1/2 (κ<0, l+1+κ=0)
        // and are present for j=l-1/2 (κ>0, l+1+κ=2l+1), splitting e.g. 2p1/2 from 2p3/2.
        double ab=es[i]+es[j];
        double t=es[i]*es[j]*::qchem::Slater::Integral(ab,2*l-1);
        if (Getκ()>0)
        {
            double kt=l+1+Getκ();
            t += kt*kt*::qchem::Slater::Integral(ab,2*l-3) - kt*ab*::qchem::Slater::Integral(ab,2*l-2);
        }
        return t*ns[i]*ns[j]; //Already has 4*Pi
    }

    virtual rvec_t     operator() (const rvec3_t&) const override;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const override;

    virtual std::string Name() const override;
private:
    rvec_t eval(const rvec3_t&) const;
};

static_assert(isOpr_Evaluator<RKBS_Evaluator>);
static_assert( is1E_Evaluator<RKBS_Evaluator>);

class Slater_Cache4 : public  Cache4
{
public:
    // using Evaluator_t = Evaluator;
    virtual void Register(Cache4_Client * eval)
    {
        assert(eval);
        auto hfeval=dynamic_cast<HF_Evaluator*>(eval);
        assert(hfeval);
        hfeval->Register(&grouper);
        //
        //  At this point we need sweep through all Cacheable4* (Rks) in Cache4::cache_t
        //  and check if geval is supported (geval.l <= Rk.LMax).
        //  All unsupport Rks will be removed.  These will then automatically be recreated next time
        //  loop_4 is called.
        //
        ::qchem::Cache4::Register(eval);
    }
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const
    {
        return new ::qchem::Slater::RkEngine(
            grouper.unique_esv[ia]+grouper.unique_esv[ib],
            grouper.unique_esv[ic]+grouper.unique_esv[id],
            grouper.LMax(ia,ib,ic,id));
    }
private:
    friend class ::Cache4Tests;
    ExponentGrouper grouper;
};


} // namespace