// File: BasisSet/Molecule/Evaluators/PG_Cart_MnD/Internal/Hermite.C
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
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Hermite;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.Types;

typedef double Array2D[LMAX+1][LMAX+1];   // Hermite1 storage (module-internal, not exported)

export namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
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
class Hermite2
{
public:
    Hermite2(double αₚ, const rvec3_t& PA, const rvec3_t& PB, int LA, int LB);
    ~Hermite2();
    double operator()(const Polarization& P,const Polarization& Pa,const Polarization& Pb) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual Hermite2* Clone(        ) const;

private:
    Hermite2(const Hermite2&) {};
    Hermite2& operator=(const Hermite2&) {return *this;}
    rvec3_t Get(int N,int na,int nb) const
    {
        size_t index=GetIndex(N,na,nb);
        return rvec3_t(d[index],e[index],f[index]);
    }
    void Assign(int N,int na,int nb,const rvec3_t& a);
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
    size_t GetSize() const {return (LA+LB)*(LA+1)*(LB+1)+LA*(LB+1)+LB+1;}

    int LA, LB;
    int LAB,LB1;
    std::vector<double> d,e,f;
    //#define USE_CACHE   // debug toggle: cache the (N,na,nb) -> index map
    #ifdef USE_CACHE
    std::map<Polarization,size_t> indexCache;
    #endif
};

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
    typedef double Array4D[3*LMAX+1][LMAX+1][LMAX+1][LMAX+1];
    friend std::ostream& operator<<(std::ostream&,const Hermite3&);

    double GetAny(const Array4D, int N, int na, int nb, int nc) const;
    double Getd(int N, int na, int nb, int nc) const {return GetAny(d,N,na,nb,nc);}
    double Gete(int N, int na, int nb, int nc) const {return GetAny(e,N,na,nb,nc);}
    double Getf(int N, int na, int nb, int nc) const {return GetAny(f,N,na,nb,nc);}

    double   itsa12s[3*LMAX+1];
    Array4D  d,e,f;
    int      itsLA, itsLB, itsLC;
    double   itsScale;
};

} //namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
