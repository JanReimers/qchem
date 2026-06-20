// File: BasisSet/DB_Cache1.C Global integrals cache allow data sharing between separate runs.
module;
#include <string>
#include <variant>
#include <functional>
export module qchem.BasisSet.Internal.DB_Cache;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI3;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Internal.Cache2;
import qchem.Types;
 
export namespace BasisSet {

// A cache client supplies its own identity string.  The cache uses it verbatim as the per-basis
// key axis and knows nothing about what is being cached (atoms, molecules, solids).  The contract:
// equal physics (same exponents / angular momenta / contraction AND same centres / orientation)
// MUST give equal strings, and any difference MUST give different strings -- the client owns ID
// assembly because it differs completely for atoms, molecules and solids.  (See IrrepBasisSet_IDs
// for the default atom assembly and PGData for the molecular, geometry-aware one.)
struct DBCacheClient
{
    virtual ~DBCacheClient() = default;
    virtual std::string BasisSetID() const = 0;
};

// Non-template bass class helps avoid so many annoying using statements in derived classes
// (in C++ typedefs don't get pulled in from template base classes)
struct IntegralsCache_Base
{
public:
    virtual ~IntegralsCache_Base() {};
    enum class I1C  {Charge, Normalization};
    enum class I2C  {Overlap, Repulsion,Kinetic,Nuclear,RestMass, InvOverlap, InvRepulsion};
    enum class I2n  {Nuclear};
    enum class I2x  {Kinetic, Repulsion,Overlap};
    enum class I3C  {Overlap, Repulsion}; // <ab|c> and <ar|1/r12|c>
    enum class I4C  {Direct,Exchange}; // <ab|cd> and <ar|1/r12|cd>
    using IBS_ID_t=std::string;        // internal per-basis key axis = DBCacheClient::BasisSetID()
    using Cluster_ID_t=std::string;
    using Mesh_ID_t=std::string;
    using RadialTypeID_t=std::string;

};

//
// Cache all integrals required for any calcuation.  Basis sets are expected to identify
// based on a RadialID and an AngularID.  Final integrals depend on both IDs.  But
// for many calculations a lot of info required for 3 and 4 center integrals depends only
// on the radial aspect of the IBSs involved.  In particular for atoms the Slater integrals
// Rk(abcd) are radial only.  Hence the separation of radial and angular IDs
//
template  <class T> class IntegralsCache : public virtual IntegralsCache_Base
{
public:

    virtual ~IntegralsCache() {};

    // Self-contained lookup-or-compute.  One call replaces the old Has()/GetXXX()/Set()
    // protocol: on a miss the make() lambda is invoked (it may itself perform nested
    // cached Get()s safely), the result is stored, and a reference to the *stored*
    // object is returned.  No shared "last key"/iterator state, so it is re-entrant.
    virtual const rvec_t&    Get(I1C,const DBCacheClient*,                          std::function<rvec_t   ()> make)=0; // Charge       -> vectors
    virtual const smat_t<T>& Get(I2C,const DBCacheClient*,                          std::function<smat_t<T>()> make)=0; // 2C sym mats
    virtual const smat_t<T>& Get(I2n,const DBCacheClient*,const Cluster_ID_t&,      std::function<smat_t<T>()> make)=0; // Nuclear
    virtual const  mat_t<T>& Get(I2x,const DBCacheClient*,const DBCacheClient*,     std::function< mat_t<T>()> make)=0; // cross IBS
    virtual const ERI3  <T>& Get(I3C,const DBCacheClient*,const DBCacheClient*,     std::function<ERI3  <T>()> make)=0; // 3 centre
    virtual const ERI4&      Get(I4C,const DBCacheClient*,const DBCacheClient*,     std::function<ERI4     ()> make)=0; // 4 centre
    // Numerically integrated variants, keyed also by a mesh ID.
    virtual const rvec_t&    Get(I1C,const DBCacheClient*,const Mesh_ID_t&,                      std::function<rvec_t()> make)=0; // Norm
    virtual const rmat_t&    Get(I2x,const DBCacheClient*,const DBCacheClient*,const Mesh_ID_t&, std::function<rmat_t()> make)=0; // mesh overlap

    // 4 center radial Slater integrals, Rk for HF calculations
    virtual void Register(Cache4_Client* eval)=0;
    virtual const Cache4* GetCache4(const RadialTypeID_t& type) const=0;

    // 2 center charge distributions Omega_ab (molecular Gaussians); analogue of Cache4.
    virtual void Register(Cache2_Client* eval)=0;
    virtual const Cache2* GetCache2(const RadialTypeID_t& type) const=0;

};

IntegralsCache<double>* theGlobalCache=0;




} //export


