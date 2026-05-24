// File: BasisSet1/Internal/DB_Cache_RAM.C  RAM implementation of database
module;
#include <map>
#include <cassert>
#include <vector>
#include <string>
#include <variant>
#include <fstream>
#include <memory>
export module qchem.BasisSet.Internal.DB_Cache_RAM;
export import qchem.BasisSet.Internal.DB_Cache;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI3;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet.Internal.Cache4;
import qchem.Types;
 
export namespace BasisSet {


template  <class T> struct IntegralsCache_RAM 
    : public virtual IntegralsCache<T>
    , public virtual IntegralsCache_Base //Get all the typedefs
{
public:
    IntegralsCache_RAM(bool makelog=false);
    virtual ~IntegralsCache_RAM(); //Report some RAM usage
    virtual bool Has(Ix1,const IBS_ID_t&) const;
    virtual bool Has(Ix2,const IBS_ID_t&,const IBS_ID_t&) const;
    virtual bool Has(I2n,const IBS_ID_t&,const Cluster_ID_t&) const;
    // Numerically integrated, with a corresponding mesh ID 
    virtual bool Has(I1C,const IBS_ID_t&,const Mesh_ID_t&   ) const;
    virtual bool Has(I2x,const IBS_ID_t&,const IBS_ID_t&,const Mesh_ID_t&) const;

    virtual const rvec_t   & GetVec () const; 
    virtual const smat_t<T>& GetSMat() const; 
    virtual const  mat_t<T>& GetMat () const; 
    virtual const ERI3  <T>& GetERI3() const; 
    virtual const ERI4     & GetERI4() const; 

    virtual const rvec_t&    Set(const rvec_t   &); 
    virtual const smat_t<T>& Set(const smat_t<T>&); 
    virtual const  mat_t<T>& Set(const  mat_t<T>&); 
    virtual const ERI3  <T>& Set(const   ERI3<T>&); 
    virtual const ERI4     & SetDirect  (const ERI4&); 
    virtual const ERI4     & SetExchange(const ERI4&); 

    void Register(Cache4_Client* eval);
    const Cache4* GetCache4(const RadialTypeID_t& type) const;

private:
    using key1_t=std::tuple<I1C,IBS_ID_t>; //Integral key for one IBS, 1 centers.
    using key2_t=std::tuple<I2C,IBS_ID_t>; //Integral key for one IBS, 2 centers.
    using keyx_t=std::tuple<I2x,IBS_ID_t,IBS_ID_t>; //Integral key for cross IBS integrals.
    using keyn_t=std::tuple<IBS_ID_t,Cluster_ID_t>; //Integral key for nuclear integrals.
    using key3_t=std::tuple<I3C,IBS_ID_t,IBS_ID_t>; //Integral key for 3 center ERI integrals between 2 IBSs.

    using key1m_t=std::tuple<I1C,IBS_ID_t,Mesh_ID_t>; //Integral key for one IBS and a mesh, 1 centers.
    using key2xm_t=std::tuple<I2x,IBS_ID_t,IBS_ID_t,Mesh_ID_t>; //Integral key for one IBS and a mesh, 1 centers.

    using map1_t=std::map<key1_t ,rvec_t   >;
    using map2_t=std::map<key2_t ,smat_t<T>>;
    using mapn_t=std::map<keyn_t ,smat_t<T>>;
    using mapx_t=std::map<keyx_t , mat_t<T>>;
    using map3_t=std::map<key3_t ,  ERI3<T>>;
    using map4_t=std::map<IBS_ID_t,std::map<IBS_ID_t,ERI4>>;

    using map1m_t =std::map<key1m_t  ,rvec_t>;
    using map2xm_t=std::map<key2xm_t ,rmat_t>;
    void ReportRAMUsage() const;
    static size_t Report(const map4_t&, const std::string&, bool verbose);

    using val_t=std::unique_ptr<Cache4>;
    using cach4_t=std::map<RadialTypeID_t,val_t>;


    mutable std::variant<typename map1_t::const_iterator,typename map1m_t::const_iterator> its1CIterator;
    mutable std::variant<typename map2_t::const_iterator,typename mapn_t::const_iterator> its2CIterator;
    mutable mapx_t::const_iterator its2xIterator;
    mutable map3_t::const_iterator its3CIterator;
    mutable std::map<IBS_ID_t,ERI4>::const_iterator its4CIterator; //Iterator into the inner map.
    mutable map2xm_t::const_iterator its2xmIterator;
    
    mutable key1_t   itsLastKey1;
    mutable key1m_t  itsLastKey1m;
    mutable key2_t   itsLastKey2;
    mutable keyn_t   itsLastKeyn;
    mutable keyx_t   itsLastKeyx;
    mutable key2xm_t itsLastKey2xm;
    mutable key3_t   itsLastKey3;
    mutable IBS_ID_t itsLastKey4a,itsLastKey4b;

    map1_t itsVecs;  //Vectors (charge)
    map2_t itsSMats; //Symmetric 2 center matrices
    mapn_t itsNMats; //Symmetric 2 center matrices for nuclear attraction integrals.
    mapx_t itsMats;  //Non-symmetric cross integrals between 2 IBSs.
    map3_t itsERI3s; //3 center, 2 IBS ERI integrals for DFT.
    map4_t Jac,Kab;  //4 center, 2 IBS ERI integrals for HF.

    map1m_t  itsmVecs; //Numerically integrated
    map2xm_t itsmMats;  //Numerically integrated

    cach4_t itsCache4s; //4 index Radial integral caches.  String identifies IBS type {Slater,BSpline<K>,POlGaussian,etc}

    bool itsMakeLog;
    mutable std::ofstream itsLogger;
};




} //namespace