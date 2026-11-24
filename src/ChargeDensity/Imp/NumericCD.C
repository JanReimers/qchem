// File: ChargeDensity/Imp/NumericCD.C  Superposition-of-atomic-densities DFT seed (impl).
module;
#include <cassert>
#include <vector>
#include <memory>
#include <cstddef>
#include <stdexcept>
#include <string>     // std::to_string (the refused count)
#include <limits>

module qchem.ChargeDensity.NumericCD;
import qchem.Blaze;                        // rvec3_t vector arithmetic (Gradient accumulation) + rvec_t
import qchem.Mesh.Quadrature;                   // qcMesh::Overlap (OverlapWith)

namespace qchem::ChargeDensity
{

NumericCD::NumericCD(double totalCharge, std::shared_ptr<const qcMesh::Mesh> quadrature)
    : itsQuadrature(std::move(quadrature)), itsCharge(totalCharge), itsVersion(NextDensityVersion())   // shared global clock (no cross-kind collisions)
{
    // An EMPTY channel is a density (a one-electron doublet's minority channel: H, Na q1); a negative or non-finite count is not, and no
    // caller can act on one (the seed computed it), so it throws -- an assert would be dead in Release (D-SEED-EMPTY-CHANNEL, 2026-11-23).
    if (!(totalCharge>=0.0) || totalCharge==std::numeric_limits<double>::infinity())
        throw std::logic_error("NumericCD: a channel's electron count must be finite and >= 0, got "+std::to_string(totalCharge));
    if (!itsQuadrature || itsQuadrature->size()==0)
        throw std::logic_error("NumericCD: a matrix-free seed needs the quadrature it integrates itself on (it is only a function on space)");
}

rvec_t NumericCD::OverlapWith(const VectorFunction<double>& g) const
{
    return qcMesh::Overlap(*itsQuadrature, g, static_cast<const ScalarFunction<double>&>(*this));
}

void NumericCD::Insert(std::shared_ptr<const ScalarFunction<double>> d)
{
    assert(d);
    itsDensities.push_back(std::move(d));
    itsVersion=NextDensityVersion();   // mutated -> a new logical density
}

double NumericCD::operator()(const rvec3_t& r) const
{
    double rho=0;
    for (const auto& d : itsDensities) rho += (*d)(r);
    return itsScale*rho;
}

rvec3_t NumericCD::Gradient(const rvec3_t& r) const
{
    rvec3_t g(0,0,0);
    for (const auto& d : itsDensities) g += d->Gradient(r);
    return itsScale*g;
}

void NumericCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}


PolarizedNumericCD::PolarizedNumericCD(std::unique_ptr<NumericCD> up, std::unique_ptr<NumericCD> down)
    : itsUp(std::move(up)), itsDn(std::move(down)), itsVersion(NextDensityVersion())
{
    if (!itsUp || !itsDn) throw std::logic_error("PolarizedNumericCD: both spin channels are required");
}

const tChargeDensity<double>* PolarizedNumericCD::GetChannel(const Spin& s) const
{
    if (s==Spin::Up)   return itsUp.get();
    if (s==Spin::Down) return itsDn.get();
    throw std::logic_error("PolarizedNumericCD::GetChannel: ask for Up or Down (None is the whole density, ChannelOf answers it)");
}

double  PolarizedNumericCD::operator()(const rvec3_t& r) const {return (*itsUp)(r)+(*itsDn)(r);}
rvec3_t PolarizedNumericCD::Gradient  (const rvec3_t& r) const {return itsUp->Gradient(r)+itsDn->Gradient(r);}
rvec_t  PolarizedNumericCD::OverlapWith(const VectorFunction<double>& g) const {return itsUp->OverlapWith(g)+itsDn->OverlapWith(g);}
double  PolarizedNumericCD::GetTotalCharge() const {return itsUp->GetTotalCharge()+itsDn->GetTotalCharge();}

void PolarizedNumericCD::ReScale(double factor)
{
    itsUp->ReScale(factor);
    itsDn->ReScale(factor);
    itsVersion=NextDensityVersion();
}

} //namespace
