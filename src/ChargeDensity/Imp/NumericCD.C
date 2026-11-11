// File: ChargeDensity/Imp/NumericCD.C  Superposition-of-atomic-densities DFT seed (impl).
module;
#include <cassert>
#include <vector>
#include <memory>
#include <cstddef>
#include <stdexcept>

module qchem.ChargeDensity.NumericCD;
import qchem.Blaze;                        // rvec3_t vector arithmetic (Gradient accumulation) + rvec_t
import qchem.Mesh.Quadrature;                   // qcMesh::Overlap (OverlapWith)

namespace qchem::ChargeDensity
{

NumericCD::NumericCD(double totalCharge, std::shared_ptr<const qcMesh::Mesh> quadrature)
    : itsQuadrature(std::move(quadrature)), itsCharge(totalCharge), itsVersion(NextDensityVersion())   // shared global clock (no cross-kind collisions)
{
    assert(totalCharge>0);
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


} //namespace
