// File: Imp/ExchangeFunctional.C   Exchange potential for DFT.
module;
#include <cassert>
module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.ChargeDensity;

namespace qchem::Hamiltonian
{

ExFunctional::ExFunctional()
    : itsChargeDensity(0)
    , isPolarized(true)
{};


void ExFunctional::InsertChargeDensity(const rChargeDensity* cd)
{
    assert(cd);
    itsChargeDensity=cd;
}

} //namespace
