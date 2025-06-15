// File: TBasisSetImplementation.C


#include "Imp/BasisSet/Fit_IBS_Common.H"
#include "Mesh/MeshIntegrator.H"
#include "oml/matrix.h"
#include "oml/vector.h"

Fit_IBS_Common::Vec Fit_IBS_Common::MakeNorm   (const Mesh* m) const
{
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Normalize(*this);
}
Fit_IBS_Common::Vec Fit_IBS_Common::MakeCharge (const Mesh*  m) const
{
    assert(false);
    return *new Vec();
}
Fit_IBS_Common::Mat Fit_IBS_Common::MakeOverlap(const Mesh* m,const fbs_t& b) const
{
    assert(false);
    return *new Mat();
}
const Fit_IBS_Common::Vec Fit_IBS_Common::Overlap  (const Mesh* m,const Sf& f) const
{
    const Vec& n=Norm(m);
    MeshIntegrator<double> mintegrator(m);
    return DirectMultiply(mintegrator.Overlap(f,*this),n);
}  
const Fit_IBS_Common::Vec Fit_IBS_Common::Repulsion(const Mesh* m,const Sf& f) const
{
    const Vec& n=Norm(m);
    MeshIntegrator<double> mintegrator(m);
    return DirectMultiply(mintegrator.Repulsion(f,*this),n);
}


