// File: BasisSet/Gaussian/Evaluators/PG_Cart_MnD/Internal/Hermite.C
//
// The Cartesian McMurchie-Davidson Hermite-coefficient blocks for the PG integrals, gathered into one
// module (the Imp stays split across Imp/Hermite1.C, Imp/Hermite2.C, Imp/Hermite3.C):
//   Hermite1 -- 1-function Hermite expansion coefficients.
//   Hermite2 -- 2-function block (a primitive pair; the Ω charge distribution's coefficients).
//   Hermite3 -- 3-function block (a primitive triple; the 3-centre overlap).
// All three are concrete value classes.  Each op() returns the product over the three Cartesian
// directions (cheap zero checks first).
module;
#include <iosfwd>
#include <vector>
#include <cassert>
#include <map>

#define LMAX 3
export module qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Internal.Hermite;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization;
import qchem.Types;

typedef double Array2D[LMAX+1][LMAX+1];   // Hermite1 storage (module-internal, not exported)

export namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{

//   n  l  m
//  d  e  f       1-function Hermite expansion coefficients; op() = product over the three directions.
//   N  L  M
class Hermite1
{
public:
    Hermite1();
    Hermite1(double αₚ, int L);

    double  operator()(const Polarization& P,const Polarization& p) const;

    void    Add  (const Hermite1&, double Scale);
    void    Clear();

private:
    double Getdef(int N,int n) const;
    int     itsL;
    Array2D def;
};

//   __  __  __
//   nn  ll  mm    2-function block (primitive pair).  GetIndex is hot in profiling.
//  d   e   f
//   N   L   M
//!
//! The constructor and \c operator() are INLINE (explicitly: a module's in-class body is not importable otherwise): the cache-free
//! 2-centre kernels in GaussianRF.C build a Hermite2 per primitive pair, and Enzyme differentiates them through it (the lattice S / T
//! pullback, Imp/GaussianRF_AD.C), so it must see these bodies.
class Hermite2
{
public:
    inline Hermite2(double αₚ, const rvec3_t& PA, const rvec3_t& PB, int LA, int LB);
    inline ~Hermite2() {}
    inline double operator()(const Polarization& P,const Polarization& Pa,const Polarization& Pb) const;

    //! Heap bytes owned by this block (the d/e/f coefficient vectors) -- for the Cache2 RAM accounting;
    //! the fixed sizeof(Hermite2) is added by the owner (Ω).
    size_t HeapBytes() const {return (d.size()+e.size()+f.size())*sizeof(double);}

    std::ostream&  Write(std::ostream&) const;   // not virtual: no class derives, and a vtable in an Enzyme-differentiated object is a hazard
    Hermite2* Clone(        ) const;

private:
    Hermite2(const Hermite2&) {};
    Hermite2& operator=(const Hermite2&) {return *this;}
    inline rvec3_t Get(int N,int na,int nb) const
    {
        size_t index=GetIndex(N,na,nb);
        return rvec3_t(d[index],e[index],f[index]);
    }
    inline void Assign(int N,int na,int nb,const rvec3_t& a)
    {
        size_t index=GetIndex(N,na,nb);
        #ifdef USE_CACHE
        Polarization key(N,na,nb);
        if (auto i=indexCache.find(key);i==indexCache.end())
            indexCache[key]=index;
        #endif
        d[index]=a.x;
        e[index]=a.y;
        f[index]=a.z;
    }
    static inline rvec3_t DirectMultiply(const rvec3_t& a, const rvec3_t& b) {return rvec3_t(a.x*b.x, a.y*b.y, a.z*b.z);}
    inline size_t GetIndex(int N, int na, int nb) const
    {
        assert(N>=0);
        assert(N<=LA+LB);
        assert(na>=0);
        assert(na<=LA);
        assert(nb>=0);
        assert(nb<=LB);
        assert(N<=na+nb);
        return N*LAB+na*LB1+nb;
    }
    inline size_t GetSize() const {return (LA+LB)*(LA+1)*(LB+1)+LA*(LB+1)+LB+1;}

    int LA, LB;
    int LAB,LB1;
    std::vector<double> d,e,f;
    //#define USE_CACHE   // debug toggle: cache the (N,na,nb) -> index map
    #ifdef USE_CACHE
    std::map<Polarization,size_t> indexCache;
    #endif
};

// The Hermite2 bodies (moved from Imp/Hermite2.C, 2026-11-13; McMurchie and Davidson, J. Comp. Phys. 26, 218 (1978); the storage
// layout is described there).
inline Hermite2::Hermite2(double αₚ, const rvec3_t& PA, const rvec3_t& PB, int _LA, int _LB)
    : LA(_LA)
    , LB(_LB)
    , LAB((LA+1)*(LB+1))
    , LB1(LB+1)
    , d(GetSize())
    , e(GetSize())
    , f(GetSize())
{
    assert(_LA>=0);
    assert(_LB>=0);
    double a12=1.0/(2*αₚ);
    Assign(0,0,0,rvec3_t(1,1,1)); //d_0^00
    //
    // Layered approach:  na+nb = layer.
    //
    for (int layer=0;layer<=LA+LB-1;layer++)
        for (int na=0;na<=layer;na++)
        {
            int nb=layer-na;
            if (na>LA || nb>LB) continue;
            for (int N=0;N<=layer+1;N++)
            {
                assert(N<=na+nb+1);
                rvec3_t t1(0,0,0),ta(0,0,0),tb(0,0,0),t3(0,0,0);
                if (N-1>=0)
                    t1=Get(N-1,na,nb)*a12;
                if (N<=na+nb)
                {
                    ta=DirectMultiply(Get(N,na,nb),PA);
                    tb=DirectMultiply(Get(N,na,nb),PB);
                }
                if (N+1<=na+nb)
                    t3=Get(N+1,na,nb)*double(N+1);

                if (na+1<=LA)
                    Assign(N,na+1,nb,t1+ta+t3);

                if (nb+1<=LB)
                    Assign(N,na,nb+1,t1+tb+t3);
            }
        }
}

inline double Hermite2::operator()(const Polarization& P,const Polarization& Pa,const Polarization& Pb) const
{
    assert(P .GetTotalL()<=LA+LB);
    assert(Pa.GetTotalL()<=LA      );
    assert(Pb.GetTotalL()<=LB      );
//
//  Take care of negative n's generated by derivatives of distributions.
//
    if ( Pa.n < 0 || Pb.n < 0 || P.n < 0) return 0.0;
    if ( Pa.l < 0 || Pb.l < 0 || P.l < 0) return 0.0;
    if ( Pa.m < 0 || Pb.m < 0 || P.m < 0) return 0.0;

    #ifdef USE_CACHE
    auto dindex=indexCache.find(Polarization(P.n,Pa.n,Pb.n));
    auto eindex=indexCache.find(Polarization(P.l,Pa.l,Pb.l));
    auto findex=indexCache.find(Polarization(P.m,Pa.m,Pb.m));
    assert(dindex!=indexCache.end());
    assert(eindex!=indexCache.end());
    assert(findex!=indexCache.end());
    return d[dindex->second]*e[eindex->second]*f[findex->second];
    #endif
    return d[GetIndex(P.n,Pa.n,Pb.n)]*e[GetIndex(P.l,Pa.l,Pb.l)]*f[GetIndex(P.m,Pa.m,Pb.m)];
}

//    _=
//   nnn          3-function block (concrete, like Hermite1/Hermite2).  No storage optimization
//  d             (on-the-fly calculation); with Scale=Eabc*(Pi/αQ)^3/2 op() returns the overlap integral.
//   0            (Was an abstract Hermite3 interface + a GaussianH3 implementation, but the Stage-1
//                single-radial collapse left only one implementation -- merged; reintroduce a base only
//                if a second 3-block type, e.g. for PG_Spherical_MnD, actually needs one.)
class Hermite3
{
public:
    Hermite3();
    Hermite3(double αₚ, const rvec3_t& PA, const rvec3_t& PB, const rvec3_t& PC, int LA, int LB, int LC, double Scale=1.0);

    double operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const;

private:
    friend std::ostream& operator<<(std::ostream&,const Hermite3&);

    // Storage is sized from (LA,LB,LC) like Hermite2's, NOT a fixed LMAX: the operator slot (C) carries a local
    // PP's r^{2n} monomials (x^6 for a GTH C4 term), and the old fixed [3*LMAX+1][LMAX+1]^3 arrays (LMAX 3) were
    // silently overrun by them in Release -- Li GTH-PADE-q3 went 9-30 mHa non-variational (OpenWork §3, 2026-11-04).
    size_t Index(int N, int na, int nb, int nc) const
    {
        assert(N >=0 && N <=itsLA+itsLB+itsLC);
        assert(na>=0 && na<=itsLA);
        assert(nb>=0 && nb<=itsLB);
        assert(nc>=0 && nc<=itsLC);
        return ((size_t(N)*(itsLA+1)+na)*(itsLB+1)+nb)*(itsLC+1)+nc;
    }
    double GetAny(const std::vector<double>&, int N, int na, int nb, int nc) const;
    double Getd(int N, int na, int nb, int nc) const {return GetAny(d,N,na,nb,nc);}
    double Gete(int N, int na, int nb, int nc) const {return GetAny(e,N,na,nb,nc);}
    double Getf(int N, int na, int nb, int nc) const {return GetAny(f,N,na,nb,nc);}

    int                 itsLA, itsLB, itsLC;
    double              itsScale;
    std::vector<double> itsa12s;
    std::vector<double> d,e,f;
};

} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
