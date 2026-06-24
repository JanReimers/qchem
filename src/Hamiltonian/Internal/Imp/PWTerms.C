// File: Hamiltonian/Internal/Imp/PWTerms.C  Plane-wave Kohn-Sham term implementations.
module;
#include <cassert>
#include <iostream>
#include <memory>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Energy;
import qchem.ChargeDensity;

namespace qchem::Hamiltonian
{

PW_External::PW_External(const cl_t& cl)
    : cStatic_HT_Imp()
    , theStructure(cl)
{
    assert(cl->GetNumAtoms()>0);
}

// Ask the basis for the configured external (pseudo)potential matrix.  The dynamic_cast is the
// sanctioned abstract->abstract move (cobs_t = Orbital_1E_IBS<dcmplx> down to the richer abstract
// DFTPotential_IBS<dcmplx>); only a basis that supports G-space DFT assembly answers it.
chmat_t PW_External::CalculateMatrix(const cobs_t* bs, const Spin&) const
{
    auto pw=dynamic_cast<const BasisSet::DFTPotential_IBS<dcmplx>*>(bs);
    assert(pw && "PW_External requires a DFTPotential_IBS (e.g. plane-wave) basis");
    return pw->MakeExternalPotential(&*theStructure);
}

void PW_External::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    te.Een=cd->DM_Contract(this);   // integral rho V_ext (the density contracts our matrix)
}

std::ostream& PW_External::Write(std::ostream& os) const
{
    return os << "    PW external (pseudo)potential with " << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

} //namespace
