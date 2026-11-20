// File: BasisSet/Radial/Evaluators/Slater/Evaluator.C
module;
#include <iosfwd>
#include <cassert>
#include <vector>
#include <string>
#include <stdexcept>
#include "forward.H"
export module qchem.BasisSet.Radial.Evaluators.Slater.IBS;
import qchem.BasisSet.Radial.Evaluators.Internal.ExponentialEvaluator;   // a base, not vocabulary: plain import like the Gaussian sibling (V1.20d)
import qchem.BasisSet.Radial.Evaluators.Internal.NR_Angular;
import qchem.BasisSet.Radial.Evaluators.Internal.RKBL_Angular;
import qchem.BasisSet.Radial.Evaluators.Slater.Internal.Integrals;
import qchem.BasisSet.Radial.Evaluators.Slater.Internal.Rk;
import qchem.BasisSet.Radial.Evaluators.Internal.Grouper;
import qchem.BasisSet.Radial.Evaluators;
import qchem.BasisSet.Internal.Cache4;
import qchem.Symmetry.Atom.Spherical;
import qchem.IntPow;
import qchem.Blaze;

export namespace qchem::BasisSet::Radial::Evaluators::Slater
{

// Slater::Radial holds all Slater-specific 1e integrals and Rk machinery,
// shared between NR (Evaluator) and RKB large component (RKBS_Evaluator).
class Radial : public ExponentialEvaluator
{
public:
    Radial(const rvec_t& _es, const sym_t& ir, size_t ltrim=0)
        : ExponentialEvaluator(_es,ir,ltrim) 
        , l(Symmetry::Atom::Getl(ir))
        {ns=norms();}
    Radial(size_t N, double emin, double emax, const sym_t& ir, size_t ltrim=0)
        : ExponentialEvaluator(exponents(N,emin,emax,ir),ir,ltrim) 
        , l(Symmetry::Atom::Getl(ir))
        {ns=norms();}
    //! \brief A COULOMB-FIT block whose functions are spherical charge distributions of their own radial power:
    //! \f$f_i(r) = n_i\, r^{2 m_i} e^{-\zeta_i r}\f$ (\a ms[i] = m_i).  An l-shell's density is \f$r^{2l} e^{-(\zeta_a+\zeta_b) r}\f$, so the fit
    //! of a heavy atom wants m = l: an s-only set (m = 0, the classic fit) builds that shape only from near-dependent exponents (Slater U:
    //! the s block x 2 costs 32 mHa, and a densified s set stalls 1.5 mHa out; lane B 2026-11-14).  \a ir must be an s (l = 0) irrep: the
    //! functions are spherical whatever their power.  With no \a ms every formula below is the classic one, to the bit.
    Radial(const rvec_t& _es, const std::vector<int>& _ms, const sym_t& ir)
        : ExponentialEvaluator(_es,ir,0)
        , l(Symmetry::Atom::Getl(ir))
        , ms(_ms)
        {
            if (l!=0) throw std::invalid_argument("Slater::Radial: a power-m fit block is spherical; its irrep must be s, got l="+std::to_string(l));
            if (ms.size()!=es.size()) throw std::invalid_argument("Slater::Radial: "+std::to_string(ms.size())+" powers for "+std::to_string(es.size())+" exponents");
            for (int m : ms) if (m<0 || m>3) throw std::invalid_argument("Slater::Radial: a fit power m="+std::to_string(m)+" outside [0,3] (the Rk tables' LMax)");
            ns=norms();
        }

    virtual std::ostream& Write   (std::ostream&) const override;
    double Overlap(size_t i,size_t j) const
    {
        if (!ms.empty()) return ::qchem::Slater::Integral(es[i]+es[j],2*ms[i]+2*ms[j])*ns[i]*ns[j];   // two power-m fit functions
        return ::qchem::Slater::Integral(es[i]+es[j],2*l)*ns[i]*ns[j]; //Already has 4*Pi and r^2 from dr.
    } 
    double Inv_r1(size_t i,size_t j) const
    {
        return ::qchem::Slater::Integral(es[i]+es[j],2*l-1)*ns[i]*ns[j]; //Already has 4*Pi
    }
    double Overlap(size_t i,size_t j, const Radial& c, size_t ic) const
    {
        const int pc = c.ms.empty() ? c.l : 2*c.ms[ic];   // the fit function's radial power
        return ::qchem::Slater::Integral(es[i]+es[j]+c.es[ic],2*l+pc)*ns[i]*ns[j]*c.ns[ic]; //Already has 4*Pi and r^2 from dr.
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

    // The Coulomb integrals see a function as a CHARGE DISTRIBUTION r^{2 m} e^{-zeta r}: R_0(m,m,m',m') (RkEngine::DirectR0).  A classic
    // block's m is its l (only l = 0 ever fits); a power-m fit block's m is per function (\c ChargeHalfPower).
    double Repulsion(size_t i,size_t j, const Radial& c, size_t ic) const
    {
        const size_t mc=c.ChargeHalfPower(ic);
        ::qchem::Slater::RkEngine cd(es[i]+es[j],c.es[ic],max(size_t(l),mc));
        return cd.DirectR0  (l,mc)*FourPi2*ns[i]*ns[j]*c.ns[ic];
    }
    double Repulsion(size_t i,size_t j) const
    {
        const size_t mi=ChargeHalfPower(i), mj=ChargeHalfPower(j);
        ::qchem::Slater::RkEngine cd(es[i],es[j],max(mi,mj));
        return cd.DirectR0  (mi,mj)*FourPi2*ns[i]*ns[j];
    }
    double Repulsion(size_t i,size_t j, const Radial& b) const
    {
        const size_t mi=ChargeHalfPower(i), mj=b.ChargeHalfPower(j);
        ::qchem::Slater::RkEngine cd(es[i],b.es[j],max(mi,mj));
        return cd.DirectR0  (mi,mj)*FourPi2*ns[i]*b.ns[j];
    }
    double Charge(size_t i) const
    {
        if (!ms.empty()) return ::qchem::Slater::Integral(es[i],2*ms[i])*ns[i];   // 4 pi int r^{2m} e^{-zeta r} r^2 dr
        return ::qchem::Slater::Integral(es[i],l)*ns[i];
    }
    double Norm(size_t i) const
    {
        if (!ms.empty()) return 1.0/sqrt(::qchem::Slater::Integral(2*es[i],4*ms[i]));
        return 1.0/sqrt(::qchem::Slater::Integral(2*es[i],2*l));
    }
    //! The half-power m of function \a i as a Coulomb charge distribution r^{2m} e^{-zeta r}: its l in a classic block, its own m in a
    //! power-m fit block.
    size_t ChargeHalfPower(size_t i) const {return ms.empty() ? size_t(l) : size_t(ms[i]);}
    //! The integral-cache key carries the powers: two fit blocks with the same exponents and different powers are different functions.
    virtual std::string RadialID() const override
    {
        std::string id=ExponentialEvaluator::RadialID();
        if (ms.empty()) return id;
        id+=" m={";
        for (int m : ms) id+=std::to_string(m)+" ";
        return id+"}";
    }

    virtual  rvec_t Norm     () const override {return ns;}

    virtual rvec_t     operator() (const rvec3_t&) const override;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const override;

    virtual std::string Name() const override;
    virtual std::string RadialType() const override;
    virtual Cache4*    MakeCache4() const override;
    using rvec11_t=rvec11_t;
    static double direct(const Cacheable4* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        // Survey fix 2026-09-08: was an unchecked dynamic_cast dereferenced on the next line --
        // a segfault, not a diagnosis, if a Cache4 from another radial family ever arrived.
        const ::qchem::Slater::RkEngine& cd = qchem::RequireEngine<::qchem::Slater::RkEngine>(c, "Slater::Evaluator::direct");
        return cd.DirectRk  (la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable4* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        // Survey fix 2026-09-08: was an unchecked dynamic_cast dereferenced on the next line --
        // a segfault, not a diagnosis, if a Cache4 from another radial family ever arrived.
        const ::qchem::Slater::RkEngine& cd = qchem::RequireEngine<::qchem::Slater::RkEngine>(c, "Slater::Evaluator::exchange");
        return cd.ExchangeRk(la,lc,Ak); // contract over k Rk*Ak, exchange version is more complicated
    }

    int l;
    std::vector<int> ms;   //!< a power-m fit block's half-powers (empty: a classic block, every function r^l e^{-zeta r})

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
    //! A power-m Coulomb-fit block (\c Radial's power-m constructor): spherical, \a ir an s irrep.
    NR_Evaluator(const rvec_t& es, const std::vector<int>& ms, const sym_t& ir)
        : Radial(es,ms,ir)
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