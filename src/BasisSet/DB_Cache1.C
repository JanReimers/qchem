// File: BasisSet/DB_Cache1.C Global integrals cache allow data sharing between separate runs.
module;
#include <map>
#include <cassert>
#include <vector>
#include <string>
#include <variant>
export module qchem.BasisSet.DB_Cache1;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI3;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.Types;
// import qchem.Orbital_HF_IBS;
// import qchem.BasisSet.Internal.Cache4;
 
export {
//
//  Abstract base for the type being cached.
//
class Cacheable1
{
public:
    virtual ~Cacheable1() {};
    virtual std::string Name() const=0;
};
//
//  Cache object based on four unsigned integer indices.  Client code using the caching
//  should derive from this class and overload all the protected functions with function forwarding.
//  Use covariant return types for the loop_4 overload. 
//  Derived class also needs to supply a Create function.
//
class Cache41
{
public:
    size_t Register(const std::string& BF_ID); //register a BasisFunction ID and return its unique index.

    virtual ~Cache41();
    void       loop_1(size_t i1) const;
    void       loop_2(size_t i2) const;
    void       loop_3(size_t i3) const;
    virtual    const Cacheable1* loop_4(size_t i4) const;
    template <class T> const T* Tloop_4(size_t i4) const
    {
        const Cacheable1* c=loop_4(i4);
        assert(c);
        const T* Tc=dynamic_cast<const T*>(c);
        assert(Tc);
        return Tc;
    } 
    
private:
    virtual const Cacheable1* Create(size_t i1,size_t i2,size_t i3,size_t i4) const=0;

    typedef std::map<size_t,const Cacheable1*> cache_4; 
    typedef std::map<size_t,cache_4> cache_3; 
    typedef std::map<size_t,cache_3> cache_2; 
    typedef std::map<size_t,cache_2> cache_t; 
    
    mutable cache_t cache;
    mutable cache_2* i1_cache;
    mutable cache_3* i2_cache;
    mutable cache_4* i3_cache;
    mutable size_t i1,i2,i3,i4; //Current indexes into unique list.

    mutable std::map<std::string,size_t> itsUniqueBFs; //unique list of basis functions and thier indexes.
};

template  <class T> struct IntegralsCache
{
public:
    static IntegralsCache* theGlobalCache;

    enum class I1C  {Charge, Normalization};
    enum class I2C  {Overlap, Repulsion,Kinetic,Nuclear,RestMass, InvOverlap, InvRepulsion};
    enum class I2n  {Nuclear};
    enum class I2x  {Kinetic};
    enum class I3C  {Overlap, Repulsion}; // <ab|c> and <ar|1/r12|c>
    enum class I4C  {Direct,Exchange}; // <ab|cd> and <ar|1/r12|cd>
    using Ix1=std::variant<I1C,I2C>; //Integrals that depend on 1 IBS id;
    using Ix2=std::variant<I2x,I3C,I4C>; //Integrals that depend on 2 IBS ids;
    using IBS_ID_t=std::tuple<std::string,std::string>; // <RadialID,AngularID> identifies and IBS
    using Cluster_ID_t=std::string;
    virtual bool Has(Ix1,const IBS_ID_t&) const=0;
    virtual bool Has(Ix2,const IBS_ID_t&,const IBS_ID_t&) const=0;
    virtual bool Has(I2n,const IBS_ID_t&,const Cluster_ID_t&) const=0;
    
    // All of these are expected to use the iterator from the latest Has() call.
    virtual const rvec_t   & GetVec () const=0; //Always real?
    virtual const smat_t<T>& GetSMat() const=0; 
    virtual const  mat_t<T>& GetMat () const=0; 
    virtual const ERI3  <T>& GetERI3() const=0; 
    virtual const ERI4     & GetERI4() const=0; 

    virtual void Set(const rvec_t   &)=0; 
    virtual void Set(const smat_t<T>&)=0; 
    virtual void Set(const  mat_t<T>&)=0; 
    virtual void Set(const   ERI3<T>&)=0; 
    virtual void Set(const   ERI4   &)=0; 
};

template  <class T> struct IntegralsCache_RAM : public virtual IntegralsCache<T>
{
public:
    using Ix1=IntegralsCache<T>::Ix1;
    using Ix2=IntegralsCache<T>::Ix2;
    using I1C=IntegralsCache<T>::I1C;
    using I2C=IntegralsCache<T>::I2C;
    using I2n=IntegralsCache<T>::I2n;
    using I2x=IntegralsCache<T>::I2x;
    using I3C=IntegralsCache<T>::I3C;
    using I4C=IntegralsCache<T>::I4C;
    using IBS_ID_t=IntegralsCache<T>::IBS_ID_t;
    using Cluster_ID_t=IntegralsCache<T>::Cluster_ID_t;

    virtual bool Has(Ix1,const IBS_ID_t&) const;
    virtual bool Has(Ix2,const IBS_ID_t&,const IBS_ID_t&) const;
    virtual bool Has(I2n,const IBS_ID_t&,const Cluster_ID_t&) const;

    virtual const rvec_t   & GetVec () const; 
    virtual const smat_t<T>& GetSMat() const; 
    virtual const  mat_t<T>& GetMat () const; 
    virtual const ERI3  <T>& GetERI3() const; 
    virtual const ERI4     & GetERI4() const; 

    virtual void Set(const rvec_t   &); 
    virtual void Set(const smat_t<T>&); 
    virtual void Set(const  mat_t<T>&); 
    virtual void Set(const   ERI3<T>&); 
    virtual void Set(const   ERI4   &); 

private:
    using key1_t=std::tuple<I1C,IBS_ID_t>; //Integral key for one IBS, 1 centers.
    using key2_t=std::tuple<I2C,IBS_ID_t>; //Integral key for one IBS, 2 centers.
    using keyx_t=std::tuple<I2x,IBS_ID_t,IBS_ID_t>; //Integral key for cross IBS integrals.
    using keyn_t=std::tuple<IBS_ID_t,Cluster_ID_t>; //Integral key for nuclear integrals.
    using key3_t=std::tuple<I3C,IBS_ID_t,IBS_ID_t>; //Integral key for 3 center ERI integrals between 2 IBSs.

    using map1_t=std::map<key1_t ,rvec_t   >;
    using map2_t=std::map<key2_t ,smat_t<T>>;
    using mapn_t=std::map<keyn_t ,smat_t<T>>;
    using mapx_t=std::map<keyx_t , mat_t<T>>;
    using map3_t=std::map<key3_t ,  ERI3<T>>;
    using map4_t=std::map<IBS_ID_t,std::map<IBS_ID_t,ERI4>>;

    mutable map1_t::const_iterator its1CIterator;
    mutable std::variant<typename map2_t::const_iterator,typename mapn_t::const_iterator> its2CnIterator;
    mutable mapx_t::const_iterator its2xIterator;
    mutable map3_t::const_iterator its3CIterator;
    mutable std::map<IBS_ID_t,ERI4>::const_iterator its4CIterator; //Iterator into the inner map.

    
    mutable key1_t itsLastKey1;
    mutable key2_t itsLastKey2;
    mutable keyn_t itsLastKeyn;
    mutable keyx_t itsLastKeyx;
    mutable key3_t itsLastKey3;
    mutable IBS_ID_t itsLastKey4a,itsLastKey4b;

    map1_t itsVecs;  //Vectors (charge)
    map2_t itsSMats; //Symmetric 2 center matrices
    mapn_t itsNMats; //Symmetric 2 center matrices for nuclear attraction integrals.
    mapx_t itsMats;  //Non-symmetric cross integrals between 2 IBSs.
    map3_t itsERI3s;  //3 center, 2 IBS ERI integrals for DFT.
    map4_t Jac,Kab; //4 center, 2 IBS ERI integrals for HF.

    std::map<std::string,Cache41*> itsCaches; //4 index Radial integral caches.  String identifies IBS type {Slater,BSpline<K>,POlGaussian,etc}


};
//
// Cache all integrals required for any calcuation.  Basis sets are expected to identify
// based on a RadialID and an AngularID.  Final integrals depend on both IDs.  But
// for many calculations a lot of info required for 3 and 4 center integrals depends only
// on the radial aspect of the IBSs involved.  In particular for atoms the Slater integrals
// Rk(abcd) are radial only.
//

template <class T> class Integrals_Overlap 
{
    using I2C=IntegralsCache<T>::I2C;
public:
    virtual smat_t<T> MakeOverlap() const=0;
    virtual std::string  RadialID() const=0;
    virtual std::string AngularID() const=0;
    
    const smat_t<T>& Overlap() const
    {
        auto cache=IntegralsCache<T>::theGlobalCache;
        assert(cache);
        if (!cache->Has(IntegralsCache<T>::I2C::Overlap,IntegralsCache<T>::IBS_ID_t(RadialID(),AngularID())))
            cache->Set(MakeOverlap());
        return cache->GetI2C();
    }
};

} //export