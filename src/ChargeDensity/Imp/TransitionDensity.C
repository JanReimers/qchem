// File: ChargeDensity/Imp/TransitionDensity.C  The factory of the AO-matrix transition density.
module;
#include <memory>
#include <vector>
module qchem.ChargeDensity.TransitionDensity;
import qchem.ChargeDensity.Internal.TransitionDensity;   // the concrete

namespace qchem::ChargeDensity
{

template <class T> std::unique_ptr<TransitionDensity<T>>
AO_TransitionDensity_Factory(std::vector<TransitionBlock<T>> blocks, std::shared_ptr<const Symmetry::SelectionRule> rule)
{
    return std::make_unique<AO_TransitionDensity<T>>(std::move(blocks), std::move(rule));
}

template std::unique_ptr<TransitionDensity<double>> AO_TransitionDensity_Factory<double>(std::vector<TransitionBlock<double>>, std::shared_ptr<const Symmetry::SelectionRule>);
template std::unique_ptr<TransitionDensity<dcmplx>> AO_TransitionDensity_Factory<dcmplx>(std::vector<TransitionBlock<dcmplx>>, std::shared_ptr<const Symmetry::SelectionRule>);

} // namespace
