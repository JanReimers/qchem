// File: BasisSet1/Atom/Evaluators/BSpline/IBS_Evaluator.C
module;
#include <bspline/Core.h>
#include <vector>
#include <memory>
#include <iosfwd>
#include <cassert>
#include <cmath>
#include <functional>
export module qchem.BasisSet.Atom.Evaluators.BSpline.IBS_r;
import qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.BasisSet.Atom.Evaluators.Concepts;
import qchem.Symmetry.Yl;
import Common.Constants;
import qchem.BasisSet.Internal.Cache4;
import Common.IntPow;

export namespace BasisSet::Atom::Evaluators::BSpline
{

template <size_t K> class BSpline_r_IBS_Evaluator : public IBS_Evaluator
{
    typedef bspline::Spline<double, K> spline_t;
public: 
    BSpline_r_IBS_Evaluator(size_t Ngrid, double rmin, double rmax, const Irrep_QNs::sym_t& ylm);
    virtual void Register(Grouper*); //Set up unique spline or exponent indexes.
    
    virtual std::ostream& Write   (std::ostream&) const;
    virtual size_t maxSpan() const {return l<=K ? K-l : 0;}  //assume no overlap for indices separated by > maxSpan

    double Overlap(size_t i,size_t j) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return BilinearForm{IdentityOperator{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Grad2(size_t i,size_t j) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        static const auto T = -Dx<2>{};
        return BilinearForm{T}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Grad2(size_t i,size_t j,const BSpline_r_IBS_Evaluator& b) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        static const auto T = -Dx<2>{};
        return BilinearForm{T}(splines[i],b.splines[j])*FourPi*ns[i]*b.ns[j]; 
    } 
    double Inv_r1(size_t i,size_t j) const //no l dependence
    {
        const spline_t &a=splines[i], &b=splines[j];
        std::function< double (double)> xm1 = [&a,&b](double r)
        {
            assert(r!=0.0);
            return a(r)*b(r)/r;
        };
        return itsGL1D->Integrate(xm1)*FourPi*ns[i]*ns[j];
    } 
    double Inv_r2(size_t i,size_t j) const //no l dependence
    {
        const spline_t &a=splines[i], &b=splines[j];
        std::function< double (double)> xm2 = [&a,&b](double r)
        {
            assert(r!=0.0);
            return a(r)*b(r)/(r*r);
        };
        return itsGL1D->Integrate(xm2)*FourPi*ns[i]*ns[j];
    } 
    double Inv_r2(size_t i,size_t j,const BSpline_r_IBS_Evaluator& _b) const //no l dependence
    {
        const spline_t &a=splines[i], &b=_b.splines[j];
        std::function< double (double)> xm2 = [&a,&b](double r)
        {
            assert(r!=0.0);
            return a(r)*b(r)/(r*r);
        };
        return itsGL1D->Integrate(xm2)*FourPi*ns[i]*_b.ns[j];
    } 

    double Charge(size_t i) const
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return LinearForm{IdentityOperator{}}(splines[i])*ns[i]*FourPi;
    }
    double Norm(size_t i) const
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return 1.0/sqrt(BilinearForm{IdentityOperator{}}(splines[i],splines[i])*FourPi);
    }
    virtual  rvec_t Norm     () const {return ns;}


    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string RadialID () const;
    virtual std::string Name    () const;
    virtual std::string RadialType() const;
    virtual Cache4*    MakeCache4() const;
    using rvec11_t=AngularIntegrals::rvec11_t;
    static double direct(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::BSpline::RkEngine<K>* cd = dynamic_cast<const ::BSpline::RkEngine<K>*>(c);
        return cd->Coulomb_Rk(la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const ::BSpline::RkEngine<K>* cd = dynamic_cast<const ::BSpline::RkEngine<K>*>(c);
        return cd->ExchangeRk(la,lc,Ak); // contract over k Rk*Ak, exchange version is more complicated
    }

    const spline_t& operator[](int index) const {return splines[index];}

protected:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
    rvec_t norms() const; //assumes es,l are already initialized

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<double> knots;
    std::vector<spline_t> splines;
    bspline::Grid<double> itsGrid;
    std::unique_ptr<GLCache1D> itsGL1D; //We have to hold a pointer, because we don't know grid early enough in the constructor.
};


static_assert(isGeneric_Evaluator<BSpline_r_IBS_Evaluator<6>>);
static_assert(is1E_Evaluator     <BSpline_r_IBS_Evaluator<6>>);
// static_assert(isFit_Evaluator    <BSpline_IBS_Evaluator<6>>);
// static_assert(isDFT_Evaluator    <BSpline_IBS_Evaluator<6>>);
static_assert(isRKBL_Evaluator   <BSpline_r_IBS_Evaluator<6>>);
static_assert(isHF_Evaluator     <BSpline_r_IBS_Evaluator<6>>);

template <size_t K> class BSpline_r_Cache4 : public  Cache4
{
public:
    BSpline_r_Cache4(const bspline::Grid<double>& grid) 
    : wp([](double r2,size_t k) {return intpow(r2,k);})
    , wm([](double r2,size_t k) {return intpow(r2,-1-k);})
    , itsMaxl(0)
    , itsGL1D(grid,K+1)
    , itsGL2D(itsGL1D,K+1)
    , itsRkCache(0) 
    {
    };
    ~BSpline_r_Cache4() {delete itsRkCache;}
    // using IBS_Evaluator_t = Gaussian_IBS_Evaluator;
    virtual void Register(Cache4_Client * eval)
    {
        assert(eval);
        BSpline_r_IBS_Evaluator<K>* geval=dynamic_cast<BSpline_r_IBS_Evaluator<K>*>(eval);
        geval->Register(&grouper);
        if (geval->Getl()>itsMaxl) itsMaxl=geval->Getl();
        //
        //  At this point we need sweep through all Cacheable* (Rks) in Cache4::cache_t
        //  and check if geval is supported (geval.l <= Rk.LMax).
        //  All unsupport Rks will be removed.  These will then automatically be recreated next time
        //  loop_4 is called.
        //
        Cache4::Register(eval);

        delete itsRkCache;
        itsRkCache=new ::BSpline::RkCache<K>(grouper.unique_spv,itsGL1D, itsMaxl,wp,wm);
    }
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const
    {
         assert(itsRkCache);
        // std::cout << "ia,ib,ic,id=" << ia << " " << ib << " " << ic << " " << id << std::endl;
        size_t lmax=grouper.LMax(ia,ib,ic,id);
        return new ::BSpline::RkEngine(grouper.unique_spv,ia,ib,ic,id,lmax,itsGL1D,itsGL2D,*itsRkCache,wp,wm);
    }
    virtual size_t RAMsize() const
    {
        size_t ndoubles=Cache4::RAMsize();
        ndoubles+=itsGL1D.RAMsize();
        ndoubles+=itsGL2D.RAMsize();
        ndoubles+=itsRkCache->RAMsize();
        return ndoubles;
    }

private:
    using func_t=::BSpline::RkEngine<K>::func_t;
    func_t wp,wm; //Weight functions for Slater integrals.
    size_t itsMaxl;
    GLCache1D   itsGL1D;
    GLCache2D   itsGL2D;

    SplineGrouper<K> grouper;
    ::BSpline::RkCache<K>* itsRkCache;
};

} //namespace