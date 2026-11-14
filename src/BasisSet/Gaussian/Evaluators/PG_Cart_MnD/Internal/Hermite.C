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
    //! Constructor and evaluator INLINE (a module's out-of-line body is not importable): the cache-free 3-centre overlap twin Enzyme
    //! differentiates (GaussianRF::ThreeCentreGradient) builds one per primitive triple.
    inline Hermite3(double αₚ, const rvec3_t& PA, const rvec3_t& PB, const rvec3_t& PC, int LA, int LB, int LC, double Scale=1.0);

    inline double operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const;

private:
    friend std::ostream& operator<<(std::ostream&,const Hermite3&);

    // Storage is sized from (LA,LB,LC) like Hermite2's, NOT a fixed LMAX: the operator slot (C) carries a local
    // PP's r^{2n} monomials (x^6 for a GTH C4 term), and the old fixed [3*LMAX+1][LMAX+1]^3 arrays (LMAX 3) were
    // silently overrun by them in Release -- Li GTH-PADE-q3 went 9-30 mHa non-variational (OpenWork §3, 2026-11-04).
    inline size_t Index(int N, int na, int nb, int nc) const
    {
        assert(N >=0 && N <=itsLA+itsLB+itsLC);
        assert(na>=0 && na<=itsLA);
        assert(nb>=0 && nb<=itsLB);
        assert(nc>=0 && nc<=itsLC);
        return ((size_t(N)*(itsLA+1)+na)*(itsLB+1)+nb)*(itsLC+1)+nc;
    }
    inline double GetAny(const std::vector<double>&, int N, int na, int nb, int nc) const;
    inline double Getd(int N, int na, int nb, int nc) const {return GetAny(d,N,na,nb,nc);}
    inline double Gete(int N, int na, int nb, int nc) const {return GetAny(e,N,na,nb,nc);}
    inline double Getf(int N, int na, int nb, int nc) const {return GetAny(f,N,na,nb,nc);}

    int                 itsLA, itsLB, itsLC;
    double              itsScale;
    std::vector<double> itsa12s;
    std::vector<double> d,e,f;
};

// Hermite3's bodies (moved from Imp/Hermite3.C 2026-11-13; the overflow prints went with them, the asserts stay).
inline Hermite3::Hermite3(double αₚ, const rvec3_t& PA, const rvec3_t& PB, const rvec3_t& PC, int LA, int LB, int LC, double Scale)
    : itsLA   (LA)
    , itsLB   (LB)
    , itsLC   (LC)
    , itsScale(Scale)
{
    const size_t n4=size_t(LA+LB+LC+1)*(LA+1)*(LB+1)*(LC+1);
    d.assign(n4,0.0); e.assign(n4,0.0); f.assign(n4,0.0);
    itsa12s.assign(LA+LB+LC+1,0.0);
    double a12=1.0/(2*αₚ);
    itsa12s[0]=1.0;
    for (int n=1; n<=itsLA+itsLB+itsLC; n++) itsa12s[n]=itsa12s[n-1]*a12;

//
//  For n<=N we can load em up.
//
    for (int n=0; n<=itsLA+itsLB+itsLC; n++)
        for (int na=0; na<=n && na<=itsLA; na++)
            for (int nb=0; nb<=n && nb<=itsLB; nb++)
                for (int nc=0; nc<=n && nc<=itsLC; nc++)
                {
                    int tn=na+nb+nc;
                    if (tn < n)
                    {
                        d[Index(n,na,nb,nc)] = 0;
                        e[Index(n,na,nb,nc)] = 0;
                        f[Index(n,na,nb,nc)] = 0;
                    }
                    if (tn == n)
                    {
                        d[Index(n,na,nb,nc)] = itsa12s[n];
                        e[Index(n,na,nb,nc)] = itsa12s[n];
                        f[Index(n,na,nb,nc)] = itsa12s[n];
                    }
                }

    for (int n=0; n<=itsLA+itsLB+itsLC; n++)
    {
//
//  For N<n we need to use the recusion relations.
//
        for (int N=0; N<n ; N++)
        {

            int na_low=n-itsLB-itsLC;
            if (na_low<0) na_low=0;
            for (int na=na_low; na<=n && na<=itsLA; na++)
            {
                int nb_low=n-na-itsLC;
                if (nb_low<0) nb_low=0;
                for (int nb=nb_low; nb<=n-na && nb<=itsLB; nb++)
                {
                    int nc=n-na-nb;


                    if (na > 0)
                    {
                        d[Index(N,na,nb,nc)] = Getd(N-1,na-1,nb,nc)*a12 + Getd(N,na-1,nb,nc)*PA.x + Getd(N+1,na-1,nb,nc)*(N+1);
                        e[Index(N,na,nb,nc)] = Gete(N-1,na-1,nb,nc)*a12 + Gete(N,na-1,nb,nc)*PA.y + Gete(N+1,na-1,nb,nc)*(N+1);
                        f[Index(N,na,nb,nc)] = Getf(N-1,na-1,nb,nc)*a12 + Getf(N,na-1,nb,nc)*PA.z + Getf(N+1,na-1,nb,nc)*(N+1);
                    }
                    else //na = 0.
                    {
                        if (nb > 0)
                        {
                            d[Index(N,na,nb,nc)] = Getd(N-1,na,nb-1,nc)*a12 + Getd(N,na,nb-1,nc)*PB.x + Getd(N+1,na,nb-1,nc)*(N+1);
                            e[Index(N,na,nb,nc)] = Gete(N-1,na,nb-1,nc)*a12 + Gete(N,na,nb-1,nc)*PB.y + Gete(N+1,na,nb-1,nc)*(N+1);
                            f[Index(N,na,nb,nc)] = Getf(N-1,na,nb-1,nc)*a12 + Getf(N,na,nb-1,nc)*PB.z + Getf(N+1,na,nb-1,nc)*(N+1);
                        }
                        else //na && nb = 0.
                        {
                            if (nc > 0)
                            {
                                d[Index(N,na,nb,nc)] = Getd(N-1,na,nb,nc-1)*a12 + Getd(N,na,nb,nc-1)*PC.x + Getd(N+1,na,nb,nc-1)*(N+1);
                                e[Index(N,na,nb,nc)] = Gete(N-1,na,nb,nc-1)*a12 + Gete(N,na,nb,nc-1)*PC.y + Gete(N+1,na,nb,nc-1)*(N+1);
                                f[Index(N,na,nb,nc)] = Getf(N-1,na,nb,nc-1)*a12 + Getf(N,na,nb,nc-1)*PC.z + Getf(N+1,na,nb,nc-1)*(N+1);
                            }
                        }
                    }


                } //nb
            } //na
        } //N
    } //n
}


inline double Hermite3::GetAny(const std::vector<double>& def, int N, int na, int nb, int nc) const
{
//  cout << "(N,n1,n2) = (" << N << "," << na << "," << nb << ")" << std::endl;
    assert(N <=itsLA+itsLB+itsLC);
    assert(na<=itsLA            );
    assert(nb<=itsLB            );
    assert(nc<=itsLC            );
//
//  Take care of negative n's generated by derivatives of distributions.
//
    if ( na < 0 || nb < 0 || nc < 0 || N < 0) return 0.0;

    int n=na+nb+nc; // we now know that na, nb and nc are >=0.
//
//  All other 0.0 elements.
//
    if (N >  n) return 0.0;
//
//  Take care of trivial elemnts.  We know n1 and n2 are >=0.
//
    if (n==N) return itsa12s[n];
//
//  Use fancy indexing to get at the numbers
//
    return def[Index(N,na,nb,nc)];
}

//#if DEBUG
inline double Hermite3::operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const
{
//  cout << Pa << Pb << Pc << " " << itsLA << " " << itsLB << " " << itsLC << std::endl;
    assert(Pa.GetTotalL()<=itsLA);
    assert(Pb.GetTotalL()<=itsLB);
    assert(Pc.GetTotalL()<=itsLC);
//  cout << "Pa=" << Pa << " Pb=" << Pb << " Pc=" << Pc << " itsScale " << itsScale << " d=" << d[0][Pa.n][Pb.n][Pc.n]
//    << " e=" << e[0][Pa.l][Pb.l][Pc.l] << " f=" << f[0][Pa.m][Pb.m][Pc.m] << std::endl;
    return itsScale * d[Index(0,Pa.n,Pb.n,Pc.n)] * e[Index(0,Pa.l,Pb.l,Pc.l)] * f[Index(0,Pa.m,Pb.m,Pc.m)];
}


} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
