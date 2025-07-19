// File: TBasisSetImplementation.C

#include <iostream>
#include <cassert>
#include <memory>
#include "Fit_IBS_Common.H"

import Mesh.Integrator;
import oml;

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
Fit_IBS_Common::Mat Fit_IBS_Common::MakeOverlap(const Mesh* m,const Fit_IBS& b) const
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


