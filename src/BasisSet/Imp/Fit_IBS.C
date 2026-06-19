// File: BasisSet/Imp/Fit_IBS.C  Imple,ent a few function for a fit basis set
module;
#include <cassert>
module qchem.BasisSet.Fit_IBS;
import qchem.Mesh.Integrator;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace BasisSet
{

const  rvec_t& Fit_IBS::Charge   () const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I1C::Charge,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),
        [this]{ return MakeCharge(); });
}
const rsmat_t& Fit_IBS::Repulsion() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::Repulsion,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),
        [this]{ return MakeRepulsion(); });
}
const  rmat_t& Fit_IBS::Repulsion(const Fit_IBS& b) const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2x::Repulsion
            ,IntegralsCache_Base::IBS_ID_t(  RadialID(),  AngularID())
            ,IntegralsCache_Base::IBS_ID_t(b.RadialID(),b.AngularID())
            ,[this,&b]{ return MakeRepulsion(b); });
}
const rsmat_t& Fit_IBS::InvOverlap() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::InvOverlap,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),
        [this]{ return MakeInvOverlap(); });
}
const rsmat_t& Fit_IBS::InvRepulsion() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::InvRepulsion,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),
        [this]{ return MakeInvRepulsion(); });
}

const rvec_t& Fit_IBS::Norm(const Mesh* m) const
{
    assert(m);
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I1C::Normalization,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),m->ID(),
        [this,m]{ return MakeNorm(m); });
}

const rmat_t& Fit_IBS::Overlap(const Mesh* m,const Fit_IBS& b) const
{
    assert(m);
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2x::Overlap,
        IntegralsCache_Base::IBS_ID_t(  RadialID(),  AngularID()),
        IntegralsCache_Base::IBS_ID_t(b.RadialID(),b.AngularID()),
        m->ID(),
        [this,m,&b]{ return MakeOverlap(m,b); });
}
rvec_t Fit_IBS::MakeNorm   (const Mesh* m) const
{
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Normalize(*this);
}

rmat_t Fit_IBS::MakeOverlap(const Mesh* m,const Fit_IBS& b) const
{
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Overlap(*this,b);
}

rvec_t Fit_IBS::Overlap  (const Mesh* m,const Sf& f) const
{
    const rvec_t& n=Norm(m);
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Overlap(f,*this) * n; //two column vectors, blaze should do component multiply (not vector dot).
}  
rvec_t Fit_IBS::Repulsion(const Mesh* m,const Sf& f) const
{
    const rvec_t& n=Norm(m);
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Repulsion(f,*this) * n; //two column vectors, blaze should do component multiply (not vector dot).
}

 rsmat_t Fit_IBS::MakeInvOverlap  () const
 {
    return blazem::inv(MakeOverlap());
 }
 rsmat_t Fit_IBS::MakeInvRepulsion() const
 {
   return blazem::inv(MakeRepulsion());
 }

} //namespace