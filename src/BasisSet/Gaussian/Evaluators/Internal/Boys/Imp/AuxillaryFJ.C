// File:  AuxillaryFJ.C  Class for calculating the auxiliary function.
module;
#include <cassert>
#include <stdexcept>
#include <string>

module qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ;


namespace qchem::BasisSet::Gaussian::Evaluators::Internal::Boys
{

AuxillaryFJ::AuxillaryFJ(int jMax) : itsjMax(jMax)
{
    // The Taylor step reads columns jmax..jmax+6 of the [121][23] table, so jmax+6 <= 22 is the real bound.
    static_assert(sizeof(theLookUp[0])/sizeof(double)==23, "Boys table width changed: re-derive thejMax");
    if (jMax<0 || jMax>thejMax)
        throw std::out_of_range("AuxillaryFJ: jmax="+std::to_string(jMax)+" exceeds the Boys F_j lookup table (jmax <= "
                               +std::to_string(thejMax)+", i.e. La+Lb+Lc+Ld <= "+std::to_string(thejMax)+")");
}

// GetFjAt is inline in the module interface (the cache-free nuclear kernel Enzyme differentiates evaluates it).

} //namespace
