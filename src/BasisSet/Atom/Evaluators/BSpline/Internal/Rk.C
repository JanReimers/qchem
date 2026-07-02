// File: src/BasisSet/Atom/Evaluators/BSpline/Internal/Rk.C  4 electron Charge distribution of BSpline orbitals. 
module;
#include <map>
#include <bspline/Core.h>
#include <functional>
export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.BasisSet.Internal.Cache4;
export import qchem.BasisSet.Atom.Evaluators.Internal.Rk;

export namespace qchem::BSpline
{

//
//  Calculate and store all off diagonal cell slater integrals.
//
template <size_t K> class RkCache
{
    typedef bspline::Spline<double,K> sp_t;
    typedef std::vector<double> dv_t;
public:
    using func_t=std::function< double (double,size_t )>;
    // \a grid selects which splines this cache covers: only pairs with BOTH splines on \a grid are
    // integrated (gl1 is built for that one grid).  The unique-spline list may span several grids now
    // that a single "BSpline<K>" Cache4 serves them all, so off-grid pairs are left empty (never queried:
    // an Rk's four splines always share a grid).
    RkCache(const std::vector<sp_t>& splines,const GLCache1D& gl1,size_t lmax, const func_t& wp, const func_t& wm, const bspline::Grid<double>& grid);
    const dv_t& plus (size_t ia,size_t ib) const {return itsMomentsPlus (ia,ib);}
    const dv_t& minus(size_t ia,size_t ib) const {return itsMomentsMinus(ia,ib);}
    size_t RAMsize() const;
private:
    smat_t<dv_t> itsMomentsPlus, itsMomentsMinus;
};

//
//  Calculate and store Rk slater integrals.
//
template <size_t K> class RkEngine  : public virtual Rk
{
    typedef bspline::Spline<double,K> sp_t;
public:
    using func_t=std::function< double (double,size_t )>;
    RkEngine(const std::vector<sp_t>& splines, size_t ia, size_t ib, size_t ic, size_t id, size_t LMax
        , const GLCache1D&,const GLCache2D&, const RkCache<K>&, const func_t& wp, const func_t& wm);
    double   DirectR0  () const; //R_0(la,la,lc,lc);
    virtual double DirectR0  (size_t la,size_t lc) const {return DirectR0  ();}
    virtual double DirectRk  (size_t la,size_t lc,const rvec11_t& Ak) const; //sum{k,A_k*R_k(la,la,lc,lc)};
    virtual double ExchangeRk(size_t la,size_t lb,const rvec11_t& Ak) const; //sum{k,A_k*R_k(la,lb,la,lb)};
    
    virtual size_t RAMsize() const;
protected:
    RkEngine(size_t _LMax) : itsLMax(_LMax), Rabcd_k(2*itsLMax+1,0.0) {};
    virtual size_t LMax() const {return itsLMax;}
    size_t itsLMax;
    rvec_t Rabcd_k;
};

} //namespace

