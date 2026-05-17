// File: BasisSet1/Atom/Evaluators/Slater/IBS_Evaluator.C
module;
#include <iosfwd>
#include <blaze/Math.h>
export module qchem.BasisSet1.Atom.Evaluators.Slater.IBS;
export import qchem.BasisSet1.Atom.Evaluators.Internal.Exponential_IBS_Evaluator;
import qchem.BasisSet1.Atom.Evaluators.Slater.Internal.Integrals; 
import qchem.BasisSet1.Atom.Evaluators.Slater.Internal.Rk; 

import qchem.BasisSet1.Internal.Cache4;

import Common.IntPow;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Yl;

export class Slater_IBS_Evaluator : public Exponential_IBS_Evaluator
{
public: 
    Slater_IBS_Evaluator(const rvec_t& es, int l, const is_t& mls) : Exponential_IBS_Evaluator(es,l,mls) {ns=norms();}
    Slater_IBS_Evaluator(const rvec_t& es, int l) : Slater_IBS_Evaluator(es,l,{}) {}
    Slater_IBS_Evaluator(const rvec_t& es, const Irrep_QNs::sym_t& ir) : Exponential_IBS_Evaluator(es,ir) {ns=norms();}
    Slater_IBS_Evaluator(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir) 
    : Exponential_IBS_Evaluator(exponents(N,emin,emax,ir),ir) {ns=norms();}
   
    Slater_IBS_Evaluator Rescale(double scale_factor) const
    {
        return Slater_IBS_Evaluator(scale_factor*es,0);
    }
    virtual std::ostream& Write   (std::ostream&) const;
 
    double Overlap(size_t i,size_t j) const
    {
        return Slater::Integral(es[i]+es[j],2*l)*ns[i]*ns[j]; //Already has 4*Pi and r^2 from dr.
    } 
    double Grad2(size_t i,size_t j) const
    {
        double ab=es[i]+es[j];
        double Term1=(l+1)*(l+1)*Slater::Integral(ab,2*l-2); //SlaterIntegral already has 4*Pi
        double Term2=-(l+1)*ab  *Slater::Integral(ab,2*l-1);
        double Term3=es[i]*es[j]*Slater::Integral(ab,2*l);
        return (Term1+Term2+Term3)*ns[i]*ns[j];
    } 
    double Grad2(size_t i,size_t j,const Slater_IBS_Evaluator& b) const
    {
        assert(l==b.l);
        double ab=es[i]+b.es[j];
        double Term1=(l+1)*(l+1)  *Slater::Integral(ab,2*l-2); //SlaterIntegral already has 4*Pi
        double Term2=-(l+1)*ab    *Slater::Integral(ab,2*l-1);
        double Term3=es[i]*b.es[j]*Slater::Integral(ab,2*l);
        return (Term1+Term2+Term3)*ns[i]*b.ns[j];
    } 
    double Inv_r1(size_t i,size_t j) const
    {
        return Slater::Integral(es[i]+es[j],2*l-1)*ns[i]*ns[j]; //Already has 4*Pi
    } 
    double Inv_r2(size_t i,size_t j) const
    {
        return Slater::Integral(es[i]+es[j],2*l-2)*ns[i]*ns[j]; //Already has 4*Pi
    } 
    double Inv_r2(size_t i,size_t j,const Slater_IBS_Evaluator& b) const
    {
        assert(l==b.l);
        return Slater::Integral(es[i]+b.es[j],2*l-2)*ns[i]*b.ns[j]; //Already has 4*Pi
    } 

    double Overlap(size_t i,size_t j, const Slater_IBS_Evaluator& c, size_t ic) const
    {
        return Slater::Integral(es[i]+es[j]+c.es[ic],2*l+c.l)*ns[i]*ns[j]*c.ns[ic]; //Already has 4*Pi and r^2 from dr.
    } 
    double Repulsion(size_t i,size_t j, const Slater_IBS_Evaluator& c, size_t ic) const
    {
        Slater::RkEngine cd(es[i]+es[j],c.es[ic],std::max(l,c.l));
        return cd.Coulomb_R0(l,c.l)*FourPi2*ns[i]*ns[j]*c.ns[ic];
    } 
    double Repulsion(size_t i,size_t j) const
    {
        Slater::RkEngine cd(es[i],es[j],l);
        return cd.Coulomb_R0(l,l)*FourPi2*ns[i]*ns[j];
    }
    double Repulsion(size_t i,size_t j, const Slater_IBS_Evaluator& b) const
    {
        Slater::RkEngine cd(es[i],b.es[j],std::max(l,b.l));
        return cd.Coulomb_R0(l,b.l)*FourPi2*ns[i]*b.ns[j];
    }
    double Charge(size_t i) const
    {
        return Slater::Integral(es[i],l)*ns[i];
    }
    double Norm(size_t i) const
    {
        return 1.0/sqrt(Slater::Integral(2*es[i],2*l));
    }

    virtual  rvec_t Norm     () const {return ns;}

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;
    virtual std::string RadialType() const;
    virtual Cache41*    MakeCache4() const;
protected:
    static rvec_t exponents(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir);
    rvec_t norms() const; //assumes es,l are already initialized

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

static_assert(isGeneric_Evaluator<Slater_IBS_Evaluator>);
static_assert(is1E_Evaluator     <Slater_IBS_Evaluator>);
static_assert(isFit_Evaluator    <Slater_IBS_Evaluator>);
static_assert(isDFT_Evaluator    <Slater_IBS_Evaluator>);
static_assert(isRKBL_Evaluator   <Slater_IBS_Evaluator>);


export class Slater_RKBS_IBS_Evaluator : public Slater_IBS_Evaluator
{
public:
    Slater_RKBS_IBS_Evaluator(const rvec_t& es, int _kappa, int l,const is_t& mls) : Slater_IBS_Evaluator(es,l,mls), kappa(_kappa) {ns=norms();}
    Slater_RKBS_IBS_Evaluator(const rvec_t& es, int _kappa, int l) : Slater_RKBS_IBS_Evaluator(es,_kappa,l,{}) {}
    Slater_RKBS_IBS_Evaluator(size_t N, double emin, double emax, int _kappa, int l): Slater_IBS_Evaluator(N,emin,emax,Irrep_QNs::sym_t(new Yl_Sym(0))), kappa(_kappa) {ns=norms();}
    rvec_t norms() const; //assumes es,l are already initialized

    double Inv_r1(size_t i,size_t j) const
    {
        return es[i]*es[j]*Slater::Integral(es[i]+es[j],2*l-1)*ns[i]*ns[j]; //Already has 4*Pi
    } 

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string Name() const;
private:
    rvec_t eval(const rvec3_t&) const;
    int kappa;
};

static_assert(isGeneric_Evaluator<Slater_RKBS_IBS_Evaluator>);
static_assert(is1E_Evaluator     <Slater_RKBS_IBS_Evaluator>);
static_assert(isFit_Evaluator    <Slater_RKBS_IBS_Evaluator>);
static_assert(isDFT_Evaluator    <Slater_RKBS_IBS_Evaluator>);
static_assert(isRKBL_Evaluator   <Slater_RKBS_IBS_Evaluator>);
