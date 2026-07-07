// File: BasisSet/DB_Cache1.C Global integrals cache allow data sharing between separate runs.
module;
#include <string>
#include <variant>
#include <functional>
export module qchem.BasisSet.Internal.DB_Cache;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI3;
import qchem.BasisSet.Internal.GMap;          // G_ERI3 (the reciprocal-space 3-centre tensor, plane-wave path)
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Internal.Cache2;
import qchem.BasisSet.Internal.Cache3;
import qchem.Types;
 
export namespace qchem::BasisSet {

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
    // The leading dimension this client expects its cached 2-centre matrices to have (= number of basis
    // functions).  The cache cross-checks it on every hit/insert, so a BasisSetID() that is not specific
    // enough (two differently-sized geometries colliding on one key) is caught right at the cache
    // boundary -- with operator/ID/dims in hand -- instead of as a Cholesky segfault three layers down.
    // Named CacheDim() (NOT size()) on purpose: size() collides with the VectorFunction/diamond
    // hierarchy.  Final overrider is the single bridge in IrrepBasisSet<T>.
    virtual size_t CacheDim() const = 0;
};

// Non-template bass class helps avoid so many annoying using statements in derived classes
// (in C++ typedefs don't get pulled in from template base classes)
struct IntegralsCache_Base
{
public:
    virtual ~IntegralsCache_Base() {};
    enum class I1C  {Charge, Normalization};
    // I2C::Kinetic = cache key for the <p^2>=<-nabla^2> building block (no 1/2; see Orbital_1E_IBS.C).
    enum class I2C  {Overlap, Repulsion,Kinetic,Nuclear,RestMass, InvOverlap, InvRepulsion};
    enum class I2n  {Nuclear};
    // I2x::Kinetic = cache key for the RKB (relativistic) L/S cross kinetic; see Orbital_DHF_IBS.C.
    enum class I2x  {Kinetic, Repulsion,Overlap};
    enum class I3C  {Overlap, Repulsion}; // <ab|c> and <ar|1/r12|c>
    enum class I4C  {Direct,Exchange}; // <ab|cd> and <ar|1/r12|cd>
    using IBS_ID_t=std::string;        // internal per-basis key axis = DBCacheClient::BasisSetID()
    using Structure_ID_t=std::string;
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
    // I2C/I2n are HERMITIAN (hmat_t): identical to symmetric for double (hmat_t<double>==smat_t<double>),
    // and correctly Hermitian (chmat_t) for the complex/plane-wave path.
    virtual const hmat_t<T>& Get(I2C,const DBCacheClient*,                          std::function<hmat_t<T>()> make)=0; // 2C herm mats
    virtual const hmat_t<T>& Get(I2n,const DBCacheClient*,const Structure_ID_t&,      std::function<hmat_t<T>()> make)=0; // Nuclear
    virtual const  mat_t<T>& Get(I2x,const DBCacheClient*,const DBCacheClient*,     std::function< mat_t<T>()> make)=0; // cross IBS
    virtual const ERI3  <T>& Get(I3C,const DBCacheClient*,const DBCacheClient*,     std::function<ERI3  <T>()> make)=0; // 3 centre
    // Reciprocal-space (plane-wave) 3-centre tensor <G_i G_j|G_c>: the G-space analogue of the ERI3 gather.
    // T-independent data (complex delta support + real kernel), but keyed/stored per cache instance like the
    // ERI3 variant; realized only on the dcmplx cache (the plane-wave path).  Overload resolves off the make
    // functor's return type (G_ERI3 vs ERI3<T>), so the two I3C Get()s never collide.
    virtual const G_ERI3&    Get(I3C,const DBCacheClient*,const DBCacheClient*,     std::function<G_ERI3   ()> make)=0; // 3 centre (G-space)
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

    // 3 center Hermite blocks (molecular DFT Overlap3C); analogue of Cache4.
    virtual void Register(Cache3_Client* eval)=0;
    virtual const Cache3* GetCache3(const RadialTypeID_t& type) const=0;

};

// Construct-on-first-use accessor for the one process-wide integrals cache (replaces the old raw
// `theGlobalCache` global that every main() had to `new` by hand and which leaked).  The instance is
// a function-static IntegralsCache_RAM<T>, so it is built on first use, destroyed at exit (its dtor
// still prints the RAM report), and its initialisation is thread-safe (C++11 magic statics).  The
// definition lives in Internal/Imp/DB_Cache.C, where IntegralsCache_RAM<T> is complete; only the
// instantiations explicitly listed there exist (currently just <double>).
template <class T> IntegralsCache<T>& theCache();




} //export


