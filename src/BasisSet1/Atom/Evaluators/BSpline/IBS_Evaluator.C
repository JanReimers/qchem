// File: BasisSet1/Atom/Evaluators/BSpline/IBS_Evaluator.C
module;
#include <bspline/Core.h>
#include <vector>
#include <memory>
#include <iosfwd>
#include <cassert>
#include <cmath>
export module qchem.BasisSet1.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet1.Atom.Evaluators.IBS;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet1.Atom.Evaluators.Internal.AngularIntegrals;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.Symmetry.Yl;
import Common.Constants;
import qchem.BasisSet1.Internal.Cache4;

//
//  This version is for phi(r) = sum(Bi(r),i)
// 
export template <size_t K> class BSpline_IBS_Evaluator : public IBS_Evaluator
{
    typedef bspline::Spline<double, K> spline_t;
public: 
    BSpline_IBS_Evaluator(size_t Ngrid, double rmin, double rmax, const Irrep_QNs::sym_t& ylm);
    BSpline_IBS_Evaluator(const BSpline_IBS_Evaluator& b) : BSpline_IBS_Evaluator(b.splines.size(),b.rmin,b.rmax,Irrep_QNs::sym_t(new Yl_Sym(0))) {};
    BSpline_IBS_Evaluator Rescale(double scale_factor) const
    {
        return BSpline_IBS_Evaluator(*this);
    }    
    virtual void Register(Grouper*); //Set up unique spline or exponent indexes.
    
    virtual std::ostream& Write   (std::ostream&) const;
    virtual size_t maxSpan() const {return l<=K ? K-l : 0;}  //assume no overlap for indices separated by > maxSpan

    double Overlap(size_t i,size_t j) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return BilinearForm{X<2>{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Grad2(size_t i,size_t j) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
        return BilinearForm{T}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Grad2(size_t i,size_t j,const BSpline_IBS_Evaluator& b) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
        return BilinearForm{T}(splines[i],b.splines[j])*FourPi*ns[i]*b.ns[j]; 
    } 
    double Inv_r1(size_t i,size_t j) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return BilinearForm{X<1>{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j];
    } 
    double Inv_r2(size_t i,size_t j) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return BilinearForm{IdentityOperator{}}(splines[i],splines[j])*FourPi*ns[i]*ns[j]; 
    } 
    double Inv_r2(size_t i,size_t j,const BSpline_IBS_Evaluator& b) const //no l dependence
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return BilinearForm{IdentityOperator{}}(splines[i],b.splines[j])*FourPi*ns[i]*b.ns[j]; 
    } 

    double Charge(size_t i) const
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return LinearForm{X<2>{}}(splines[i])*ns[i]*FourPi;
    }
    double Norm(size_t i) const
    {
        using namespace bspline::integration;
        using namespace bspline::operators; 
        return 1.0/sqrt(BilinearForm{X<2>{}}(splines[i],splines[i])*FourPi);
    }
    virtual  rvec_t Norm     () const {return ns;}

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string RadialID() const;
    virtual std::string Name    () const;
    virtual std::string RadialType() const;
    virtual Cache41*    MakeCache4() const;
    using rvec11_t=AngularIntegrals::rvec11_t;
    static double direct(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const BSpline::RkEngine<K>* cd = dynamic_cast<const BSpline::RkEngine<K>*>(c);
        return cd->Coulomb_Rk(la,lc,Ak); // contract over k Rk*Ak
    }
    static double exchange(const Cacheable* c, size_t la, size_t lc,const rvec11_t& Ak)
    {
        const BSpline::RkEngine<K>* cd = dynamic_cast<const BSpline::RkEngine<K>*>(c);
        return cd->ExchangeRk(la,lc,Ak); // contract over k Rk*Ak, exchange version is more complicated
    }
    const spline_t& operator[](int index) const {return splines[index];}

protected:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
    rvec_t norms() const; //assumes es,l are already initialized

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<double> knots;
    std::vector<spline_t> splines;
    std::unique_ptr<GLCache> itsGL;


};

static_assert(isGeneric_Evaluator<BSpline_IBS_Evaluator<6>>);
static_assert(is1E_Evaluator     <BSpline_IBS_Evaluator<6>>);
// static_assert(isFit_Evaluator    <BSpline_IBS_Evaluator<6>>);
// static_assert(isDFT_Evaluator    <BSpline_IBS_Evaluator<6>>);
static_assert(isRKBL_Evaluator   <BSpline_IBS_Evaluator<6>>);
static_assert(isHF_Evaluator     <BSpline_IBS_Evaluator<6>>);

export template <size_t K> class BSpline_Cache4 : public  Cache41
{
public:
    BSpline_Cache4() : itsRkCache(0) 
    {
    };
    ~BSpline_Cache4() {delete itsRkCache;}
    // using IBS_Evaluator_t = Gaussian_IBS_Evaluator;
    virtual void Register(Cache4_Client * eval)
    {
        assert(eval);
        BSpline_IBS_Evaluator<K>* geval=dynamic_cast<BSpline_IBS_Evaluator<K>*>(eval);
        geval->Register(&grouper);
        //
        //  At this point we need sweep through all Cacheable* (Rks) in Cache41::cache_t
        //  and check if geval is supported (geval.l <= Rk.LMax).
        //  All unsupport Rks will be removed.  These will then automatically be recreated next time
        //  loop_4 is called.
        //
        Cache41::Register(eval);

        delete itsRkCache;
        size_t lmax=grouper.maxl();
        itsRkCache=new BSpline::RkCache<K>(grouper.unique_spv,*grouper.GetGL(lmax),lmax);
    }
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const
    {
         assert(itsRkCache);
        // std::cout << "ia,ib,ic,id=" << ia << " " << ib << " " << ic << " " << id << std::endl;
        size_t lmax=grouper.LMax(ia,ib,ic,id);
        const GLCache* gl=grouper.GetGL(lmax);
        return new BSpline::RkEngine(grouper.unique_spv,ia,ib,ic,id,lmax,*gl,*itsRkCache);
    }
private:
   
    SplineGrouper<K> grouper;
    BSpline::RkCache<K>* itsRkCache;
};


export template <size_t K> class BSpline_r_IBS_Evaluator : public IBS_Evaluator
{
    typedef bspline::Spline<double, K> spline_t;
public: 
    BSpline_r_IBS_Evaluator(size_t Ngrid, double rmin, double rmax, const Irrep_QNs::sym_t& ylm);
    virtual void Register(Grouper*); //Set up unique spline or exponent indexes.
    
    virtual std::ostream& Write   (std::ostream&) const;
    virtual size_t maxSpan() const {return l<=K ? K-l : 0;}  //assume no overlap for indices separated by > maxSpan

    virtual rsmat_t Overlap  () const;
    virtual rsmat_t Grad2    () const;
    virtual rsmat_t Inv_r1   () const;
    virtual rsmat_t Inv_r2   () const;
    virtual rsmat_t Repulsion() const;
    virtual  rvec_t Charge   () const;
    virtual  rvec_t Norm     () const {return ns;}

    virtual dERI3  Overlap  (const IBS_Evaluator&) const; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator&) const; //3 center

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string RadialID () const;
    virtual std::string Name    () const;

    const spline_t& operator[](int index) const {return splines[index];}

protected:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
    rvec_t norms() const; //assumes es,l are already initialized

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<double> knots;
    std::vector<spline_t> splines;
    std::unique_ptr<GLCache> itsGL;
};