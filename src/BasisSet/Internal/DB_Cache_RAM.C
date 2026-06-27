// File: BasisSet/Internal/DB_Cache_RAM.C  RAM implementation of database
module;
#include <map>
#include <cassert>
#include <vector>
#include <string>
#include <variant>
#include <fstream>
#include <memory>
#include <chrono>
#include <functional>
export module qchem.BasisSet.Internal.DB_Cache_RAM;
export import qchem.BasisSet.Internal.DB_Cache;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI3;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Internal.Cache2;
import qchem.BasisSet.Internal.Cache3;
import qchem.Types;
 
export namespace BasisSet {


template  <class T> struct IntegralsCache_RAM 
    : public virtual IntegralsCache<T>
    , public virtual IntegralsCache_Base //Get all the typedefs
{
public:
    IntegralsCache_RAM(bool makelog=false);
    virtual ~IntegralsCache_RAM(); //Report some RAM usage

    virtual const rvec_t&    Get(I1C,const DBCacheClient*,                          std::function<rvec_t   ()>);
    virtual const hmat_t<T>& Get(I2C,const DBCacheClient*,                          std::function<hmat_t<T>()>);
    virtual const hmat_t<T>& Get(I2n,const DBCacheClient*,const Structure_ID_t&,      std::function<hmat_t<T>()>);
    virtual const  mat_t<T>& Get(I2x,const DBCacheClient*,const DBCacheClient*,     std::function< mat_t<T>()>);
    virtual const ERI3  <T>& Get(I3C,const DBCacheClient*,const DBCacheClient*,     std::function<ERI3  <T>()>);
    virtual const ERI4     & Get(I4C,const DBCacheClient*,const DBCacheClient*,     std::function<ERI4     ()>);
    virtual const rvec_t&    Get(I1C,const DBCacheClient*,const Mesh_ID_t&,                      std::function<rvec_t()>);
    virtual const rmat_t&    Get(I2x,const DBCacheClient*,const DBCacheClient*,const Mesh_ID_t&, std::function<rmat_t()>);

    void Register(Cache4_Client* eval);
    const Cache4* GetCache4(const RadialTypeID_t& type) const;

    void Register(Cache2_Client* eval);
    const Cache2* GetCache2(const RadialTypeID_t& type) const;

    void Register(Cache3_Client* eval);
    const Cache3* GetCache3(const RadialTypeID_t& type) const;

private:
    using key1_t=std::tuple<I1C,IBS_ID_t>; //Integral key for one IBS, 1 centers.
    using key2_t=std::tuple<I2C,IBS_ID_t>; //Integral key for one IBS, 2 centers.
    using keyx_t=std::tuple<I2x,IBS_ID_t,IBS_ID_t>; //Integral key for cross IBS integrals.
    using keyn_t=std::tuple<IBS_ID_t,Structure_ID_t>; //Integral key for nuclear integrals.
    using key3_t=std::tuple<I3C,IBS_ID_t,IBS_ID_t>; //Integral key for 3 center ERI integrals between 2 IBSs.

    using key1m_t=std::tuple<I1C,IBS_ID_t,Mesh_ID_t>; //Integral key for one IBS and a mesh, 1 centers.
    using key2xm_t=std::tuple<I2x,IBS_ID_t,IBS_ID_t,Mesh_ID_t>; //Integral key for one IBS and a mesh, 1 centers.

    using map1_t=std::map<key1_t ,rvec_t   >;
    using map2_t=std::map<key2_t ,hmat_t<T>>;
    using mapn_t=std::map<keyn_t ,hmat_t<T>>;
    using mapx_t=std::map<keyx_t , mat_t<T>>;
    using map3_t=std::map<key3_t ,  ERI3<T>>;
    using map4_t=std::map<IBS_ID_t,std::map<IBS_ID_t,ERI4>>;

    using map1m_t =std::map<key1m_t  ,rvec_t>;
    using map2xm_t=std::map<key2xm_t ,rmat_t>;

    // using time_t=std::chrono::sys_seconds;
    using time_t=decltype(std::chrono::system_clock::now());
    using id_pair_t=std::pair<IBS_ID_t,IBS_ID_t>;
    mutable std::map<id_pair_t,time_t> ERI4_timestamps; //Get functions need to update this
    // std::map<std::chrono::sys_seconds,std::pair<IBS_ID_t,IBS_ID_t>> ERI4_Oldest;  //Sort by time stamp

    void ReportRAMUsage(std::ostream&) const;
    static size_t Report(const map4_t&, const std::string&, bool verbose); // in MB
    // GC evicts the oldest ERI4s; `protect` is the entry just inserted (never evicted).
    void RunGarbageCollector(const id_pair_t& protect);
    size_t Purge(map4_t& eri4s,const id_pair_t& old,const id_pair_t& protect);

    using val_t=std::unique_ptr<Cache4>;
    using cach4_t=std::map<RadialTypeID_t,val_t>;
    using val2_t=std::unique_ptr<Cache2>;
    using cach2_t=std::map<RadialTypeID_t,val2_t>;
    using val3_t=std::unique_ptr<Cache3>;
    using cach3_t=std::map<RadialTypeID_t,val3_t>;

    map1_t itsVecs;  //Vectors (charge)
    map2_t itsSMats; //Symmetric 2 center matrices
    mapn_t itsNMats; //Symmetric 2 center matrices for nuclear attraction integrals.
    mapx_t itsMats;  //Non-symmetric cross integrals between 2 IBSs.
    map3_t itsERI3s; //3 center, 2 IBS ERI integrals for DFT.
    map4_t Jac,Kab;  //4 center, 2 IBS ERI integrals for HF.

    map1m_t  itsmVecs; //Numerically integrated
    map2xm_t itsmMats;  //Numerically integrated

    cach4_t itsCache4s; //4 index Radial integral caches.  String identifies IBS type {Slater,BSpline<K>,POlGaussian,etc}
    cach2_t itsCache2s; //2 index charge-distribution caches (molecular Omega_ab), keyed by RadialType
    cach3_t itsCache3s; //3 index Hermite-block caches (molecular DFT Overlap3C), keyed by RadialType

    bool itsMakeLog;
    mutable std::ofstream itsLogger;
    mutable size_t itsMaxRAM,itsTotalRAM; // in MB
};




} //namespace