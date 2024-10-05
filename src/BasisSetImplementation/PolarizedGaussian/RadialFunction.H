// File: RadialFunction.H  Abstract interface for the radial part of a basis function.
#ifndef _RadialFunction_H_
#define _RadialFunction_H_

#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianCD.H"
#include "Misc/Persistent/Persistent.H"
#include "Functions/ScalarFunction.H"
#include "Misc/Polarization.H"
#include "oml/vector3d.h"
#include <vector>
#include <map>

class   Polarization;
class   Hermite1;
class   Hermite3;
class   RNLM;
class  BasisFunctionBlock;
struct BasisFunctionBlockQuartet;
struct BasisFunctionBlockTriplet;
struct BasisFunctionBlockPair;

class CDcache
{
public:
    CDcache();
    ~CDcache();
    const GaussianCD& find(const GaussianRF*,const GaussianRF*);
    size_t size() const {return cache.size();}
    void Report(std::ostream&) const;
private:
    typedef std::pair<UniqueID::IDtype,UniqueID::IDtype> ids_t;
    static ids_t Sort(UniqueID::IDtype,UniqueID::IDtype);
    
    size_t lookups;
    std::map<ids_t,const GaussianCD*> cache;
};
//
//  Radial functions will be shared by many basis functions, for example 2Px, 2Py and
//  2Pz all have the same radial part, just different polarizations.
//  It is assumed that these radial functions (pointers) will stored in an optr_vector1<>
//  and the SortingIndex is the offset into that list.
//
class RadialFunction
    : public virtual Persistent
    , public virtual ScalarFunction<double>
{
public:

    enum Types3C {Overlap3C, Repulsion3C};
    enum Types2C {Overlap2C, Repulsion2C,Kinetic,Nuclear};

    typedef BasisFunctionBlock BFB;
    typedef RadialFunction     RF;
    typedef BasisFunctionBlockPair    BFBP;
    typedef BasisFunctionBlockTriplet BFBT;
    typedef BasisFunctionBlockQuartet BFBQ;

    virtual       bool      operator==      (const RadialFunction&) const=0; //Ignores L!
    virtual const RVec3&    GetCenter       (                     ) const=0;
    virtual       int       GetL            (                     ) const=0; //N+L+M?
    virtual       double    GetNormalization(const Polarization&  ) const=0;
    virtual       double    GetCharge       (const Polarization&  ) const=0;

    virtual void Get2CenterIntegrals(Types2C, BFBP&, SMat&, const Cluster*, double scale) const=0;
    virtual void Get2CenterIntegrals(Types2C, BFBP&,  Mat&,                 double scale) const=0; //Cross overlap & repulsion.
    virtual void Get3CenterIntegrals(Types3C, BFBT&, std::vector<SMat >&, double scale)=0;
    
    virtual double Integrate(Types2C,const RadialFunction* rb, const Polarization& pa, const Polarization& pb,CDcache& cache,const Cluster* cl=0) const=0;

    virtual void GetRepulsion4C(BFBQ&, ERIList&, double scale)=0;

    virtual       Hermite3* GetH3 (const RadialFunction&, const RadialFunction&) const=0;
    virtual const Hermite1& GetH1 (                                            ) const=0;
    virtual Matrix<double> GetAux(const std::vector<Polarization>& N,
                                  const std::vector<Polarization>& Pc,
                                  int LP, double AlphaP,const RVec3& P         ) const=0;

    virtual RadialFunction* Clone      (const RVec3&) const=0;
    virtual RadialFunction* Clone      (            ) const=0;
    static  RadialFunction* Factory    (std::istream&    )        ;
};

#endif //_RadialFunction_H_
