// File: VxcPol.C  Polarized HF exchange potential
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Types;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.Symmetry.Spin;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real exchange potential,  Vxc(ro(r)), where ro is the charge density.
//
// The Hamiltonain matrix elements are calculated
//             /
//  Vxc(i,j) = | dr Vxcfit(ro(r)) Oi(r) Oj(r) .
//             /
//
//           = Sum  { Ck <Oi|Vk|Oj> } .
//
// Polarized HF exchange as two spin-channel Vxc(-1): F^sigma += -K[D^sigma].  GetMatrix/GetEnergy just
// dispatch to the right sub-term, handing it THIS spin's density; each sub-Vxc does the whole-system
// contraction + caching.  (Mirrors FittedVxcPol.)  A spin-resolved density is required (the SCF density
// always resolves its channels once orbitals exist -- V1.37: through the face, never a container type).
VxcPol::VxcPol()  : itsUpVxc(new Vxc(-1.0)), itsDownVxc(new Vxc(-1.0)) {}
VxcPol::~VxcPol() { delete itsUpVxc; delete itsDownVxc; }

const rsmat_t& VxcPol::GetMatrix(const robs_t* bs,const Spin& s,const rChargeDensity* cd,const rbs_t* wholeBasis) const
{
    // A polarized term has no Spin::None block to hand back -- and the caller has to be able to SEE that
    // (R2.5: exit(-1) killed the pybind GUI and the test runner outright, taking the diagnostic with it).
    if (s==Spin::None)
        throw std::runtime_error("VxcPol::GetMatrix: asked for the Spin::None (unpolarized) block of a "
                                 "polarized exchange term -- a polarized term has an Up and a Down block, "
                                 "no total.");
    const rChargeDensity* SpinCD = ChannelOf(cd,s);   // this spin's density
    assert(SpinCD && "VxcPol: density must be polarized");
    return (s==Spin::Up ? itsUpVxc : itsDownVxc)->GetMatrix(bs,s,SpinCD,wholeBasis);
}
void VxcPol::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    // Sum K^alpha and K^beta: each sub-Vxc contracts ITS spin's density (0.5 Tr(D^sigma . -K^sigma)).
    const rDM_CD* up=DM_ChannelOf(cd,Spin::Up  );
    const rDM_CD* dn=DM_ChannelOf(cd,Spin::Down);
    assert(up && dn && "VxcPol energy: density must be polarized (D-backed channels)");
    itsUpVxc  ->GetEnergy(te, up);
    itsDownVxc->GetEnergy(te, dn);
}

std::ostream& VxcPol::Write(std::ostream& os) const
{
    os << "    Polarized Hartee-Fock exchange potential." << std::endl;
    return os;
}

} //namespace
