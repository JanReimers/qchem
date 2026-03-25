// File: Imp/Fit_IBS_Common.C Common implementation for Irrep fit basis set.
module;
#include <iostream>
#include <cassert>
#include <memory>
#include "blaze/Math.h"
module qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.Mesh.Integrator;

rvec_t Fit_IBS_Common::MakeNorm   (const Mesh* m) const
{
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Normalize(*this);
}
rvec_t Fit_IBS_Common::MakeCharge (const Mesh*  m) const
{
    assert(false);
    return *new rvec_t();
}
rmat_t Fit_IBS_Common::MakeOverlap(const Mesh* m,const Fit_IBS& b) const
{
    assert(false);
    return *new rmat_t();
}
const rvec_t Fit_IBS_Common::Overlap  (const Mesh* m,const Sf& f) const
{
    const rvec_t& n=Norm(m);
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Overlap(f,*this) * n; //two column vectors, blaze should do component multiply (not vector dot).
}  
const rvec_t Fit_IBS_Common::Repulsion(const Mesh* m,const Sf& f) const
{
    const rvec_t& n=Norm(m);
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Repulsion(f,*this) * n; //two column vectors, blaze should do component multiply (not vector dot).
}


