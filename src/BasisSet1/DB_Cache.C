// File: BasisSet/DB_Cache1.C Global integrals cache allow data sharing between separate runs.
module;
#include <string>
#include <variant>
export module qchem.BasisSet.DB_Cache;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI3;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet.Internal.Cache4;
import qchem.Types;
 
export namespace BasisSet {
    

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
    using Ix1=std::variant<I1C,I2C>; //Integrals that depend on 1 IBS id;
    using Ix2=std::variant<I2x,I3C,I4C>; //Integrals that depend on 2 IBS ids;
    using IBS_ID_t=std::tuple<std::string,std::string>; // <RadialID,AngularID> identifies and IBS
    using Cluster_ID_t=std::string;
    using Mesh_ID_t=std::string;
    using RadialTypeID_t=std::string;

    

    virtual bool Has(Ix1,const IBS_ID_t&                    ) const=0;
    virtual bool Has(Ix2,const IBS_ID_t&,const IBS_ID_t&    ) const=0;
    virtual bool Has(I2n,const IBS_ID_t&,const Cluster_ID_t&) const=0;

    // Numerically integrated, with a corresponding mesh ID 
    virtual bool Has(I1C,const IBS_ID_t&,const Mesh_ID_t&   ) const=0;
    virtual bool Has(I2x,const IBS_ID_t&,const IBS_ID_t&,const Mesh_ID_t&) const=0;
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

    // All of these are expected to use the iterator from the latest Has() or Set() call.
    virtual const rvec_t   & GetVec () const=0; //Always real?
    virtual const smat_t<T>& GetSMat() const=0; 
    virtual const  mat_t<T>& GetMat () const=0; 
    virtual const ERI3  <T>& GetERI3() const=0; 
    virtual const ERI4     & GetERI4() const=0; 

    // Insert and cache the iterator.  Return ref to the *stored* object.
    virtual const rvec_t&    Set(const rvec_t   &)=0; 
    virtual const smat_t<T>& Set(const smat_t<T>&)=0; 
    virtual const  mat_t<T>& Set(const  mat_t<T>&)=0; 
    virtual const ERI3  <T>& Set(const   ERI3<T>&)=0; 
    virtual const ERI4&      SetDirect  (const ERI4&)=0; 
    virtual const ERI4&      SetExchange(const ERI4&)=0; 

    // 4 center radial Slater integrals, Rk for HF calculations
    virtual void Register(Cache4_Client* eval)=0;
    virtual const Cache4* GetCache4(const RadialTypeID_t& type) const=0;

};

IntegralsCache<double>* theGlobalCache=0;




} //export


