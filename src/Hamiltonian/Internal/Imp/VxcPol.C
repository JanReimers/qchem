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
// contraction + caching.  (Mirrors FittedVxcPol.)  A polarized density is required (the SCF density is
// always a Polarized_CD once orbitals exist).
VxcPol::VxcPol()  : itsUpVxc(new Vxc(-1.0)), itsDownVxc(new Vxc(-1.0)) {}
VxcPol::~VxcPol() { delete itsUpVxc; delete itsDownVxc; }

const rsmat_t& VxcPol::GetMatrix(const robs_t* bs,const Spin& s,const rChargeDensity* cd,const rbs_t* wholeBasis) const
{
    if (s==Spin::None) { std::cerr << "VxcPol::GetMatrix: unpolarized spin in a polarized term" << std::endl; exit(-1); }
    const Polarized_CD* pcd = dynamic_cast<const Polarized_CD*>(cd);
    assert(pcd && "VxcPol: density must be polarized");
    const rDM_CD* SpinCD = pcd->GetChargeDensity(s);   // this spin's density
    return (s==Spin::Up ? itsUpVxc : itsDownVxc)->GetMatrix(bs,s,SpinCD,wholeBasis);
}
void VxcPol::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    // Sum K^alpha and K^beta: each sub-Vxc contracts ITS spin's density (0.5 Tr(D^sigma . -K^sigma)).
    const Polarized_CD* pcd=dynamic_cast<const Polarized_CD*>(cd);
    assert(pcd && "VxcPol energy: density must be polarized");
    itsUpVxc  ->GetEnergy(te, pcd->GetChargeDensity(Spin::Up  ));
    itsDownVxc->GetEnergy(te, pcd->GetChargeDensity(Spin::Down));
}

std::ostream& VxcPol::Write(std::ostream& os) const
{
    os << "    Polarized Hartee-Fock exchange potential." << std::endl;
    return os;
}

} //namespace
