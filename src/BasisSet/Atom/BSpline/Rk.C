// File: BasisSet/Atom/radial/BSpline/Rk.C  4 electron Charge distribution of BSpline orbitals. 
module;
#include <map>
#include <bspline/Core.h>
export module qchem.BasisSet.Atom.BSpline.Rk;
import qchem.Basisset.Atom.BSpline.GLQuadrature;
import qchem.BasisSet.Internal.Cache4;
export import qchem.BasisSet.Atom.Rk;
// import oml.Matrix; 

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
    RkCache(const std::vector<sp_t>& splines,const GLCache& gl,size_t lmax);
    const dv_t& find_plus (size_t ia,size_t ib) const {return find(ia,ib,itsMomentsPlus);}
    const dv_t& find_minus(size_t ia,size_t ib) const {return find(ia,ib,itsMomentsMinus);}
private:
    typedef std::pair<size_t,size_t> id2_t; //convention id_1 < id_2
    typedef std::map<id2_t,dv_t> moment_t;
    static const dv_t& find(size_t ia,size_t ib,const moment_t&);


    moment_t itsMomentsPlus;  //<r^(2+k) B1*B2>
    moment_t itsMomentsMinus;  //<r^(1-k) B1*B2>
};
//
// Store Slater repulsions integral for 0 <= k <= 2*LMax.
//
template <size_t K> class RkEngine  : public virtual Rk
{
    typedef bspline::Spline<double,K> sp_t;
public:
    RkEngine(const std::vector<sp_t>& splines, size_t ia, size_t ib, size_t ic, size_t id, size_t LMax, const GLCache& gl, const RkCache<K>&);
    double Coulomb_R0(size_t la,size_t lc) const {return Coulomb_R0();}
    double Coulomb_R0() const; //R_0(la,la,lc,lc);
    RVec   Coulomb_Rk(size_t la,size_t lc) const; 
    RVec   ExchangeRk(size_t la,size_t lb) const; 
private:
    size_t LMax;
    Vector<double> Rabcd_k;
};



} //namespace

