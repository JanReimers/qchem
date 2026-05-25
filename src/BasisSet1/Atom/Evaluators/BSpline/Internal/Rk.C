// File: src/BasisSet1/Atom/Evaluators/BSpline/Internal/Rk.C  4 electron Charge distribution of BSpline orbitals. 
module;
#include <map>
#include <bspline/Core.h>
export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet.Internal.Cache4;
export import qchem.BasisSet.Atom.Evaluators.Internal.Rk;

export namespace BSpline
{

//
//  Calculate and store all off diagonal and diagonal cell slater integrals.
//
template <size_t K> class RkCache
{
    typedef bspline::Spline<double,K> sp_t;
    typedef std::vector<double> dv_t;
public:
    RkCache(const std::vector<sp_t>& splines,const GLCache1D& gl1,size_t lmax);
    const dv_t& find_plus (size_t ia,size_t ib) const {return find(ia,ib,itsMomentsPlus);}
    const dv_t& find_minus(size_t ia,size_t ib) const {return find(ia,ib,itsMomentsMinus);}
    size_t RAMsize() const;
private:
    typedef std::pair<size_t,size_t> id2_t; //convention id_1 < id_2
    typedef std::map<id2_t,dv_t> moment_t;
    static const dv_t& find(size_t ia,size_t ib,const moment_t&);


    moment_t itsMomentsPlus;  //<r^(k+2) B1*B2>
    moment_t itsMomentsMinus;  //<r^(1-k) B1*B2>
};
template <size_t K> class RkCache_r
{
    typedef bspline::Spline<double,K> sp_t;
    typedef std::vector<double> dv_t;
public:
    RkCache_r(const std::vector<sp_t>& splines,const GLCache1D& gl,size_t lmax);
    const dv_t& find_plus (size_t ia,size_t ib) const {return find(ia,ib,itsMomentsPlus);}
    const dv_t& find_minus(size_t ia,size_t ib) const {return find(ia,ib,itsMomentsMinus);}
private:
    typedef std::pair<size_t,size_t> id2_t; //convention id_1 < id_2
    typedef std::map<id2_t,dv_t> moment_t;
    static const dv_t& find(size_t ia,size_t ib,const moment_t&);


    moment_t itsMomentsPlus;  //<r^k B1*B2>
    moment_t itsMomentsMinus;  //<r^(-1-k) B1*B2>
};
//
// Store Slater repulsions integral for 0 <= k <= 2*LMax.
//  This version is for phi(r) = sum(Bi(r),i)
// 
template <size_t K> class RkEngine  : public virtual Rk
{
    typedef bspline::Spline<double,K> sp_t;
public:
    RkEngine(const std::vector<sp_t>& splines, size_t ia, size_t ib, size_t ic, size_t id, size_t LMax
        , const GLCache1D&,const GLCache2D&, const RkCache<K>&);
    double   Coulomb_R0() const; //R_0(la,la,lc,lc);
    virtual double Coulomb_R0(size_t la,size_t lc) const {return Coulomb_R0();}
    virtual double Coulomb_Rk(size_t la,size_t lc,const rvec11_t& Ak) const; //sum{k,A_k*R_k(la,la,lc,lc)};
    virtual double ExchangeRk(size_t la,size_t lb,const rvec11_t& Ak) const; //sum{k,A_k*R_k(la,lb,la,lb)};
    
    virtual bool   isSupported(const Cache4_Client*) const;
    virtual size_t RAMsize() const;
protected:
    RkEngine(size_t _LMax) : LMax(_LMax), Rabcd_k(2*LMax+1,0.0) {};
    size_t LMax;
    rvec_t Rabcd_k;
};

//  This version is for phi(r) = 1/r * sum(Bi(r),i)
template <size_t K> class RkEngine_r  : public virtual Rk, private RkEngine<K>
{
    typedef bspline::Spline<double,K> sp_t;
    using RkEngine<K>::LMax;
    using RkEngine<K>::Rabcd_k;
public:
    RkEngine_r(const std::vector<sp_t>& splines, size_t ia, size_t ib, size_t ic, size_t id, size_t LMax, const GLCache1D&, const GLCache2D&, const RkCache_r<K>&);
};


} //namespace

