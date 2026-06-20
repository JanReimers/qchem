// File: DiracKinetic.C  Kinetic energy term for the Dirac hamiltonian.
module;
#include <iostream>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.Constants;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

// The RELATIVISTIC kinetic ENERGY term (Dirac \f$c\,\vec\sigma\cdot\vec p\f$, L/S off-diagonal).
// Here bs is an RKB basis, so bs->Kinetic() is NOT the bare \f$\langle p^2\rangle\f$ block -- it is the
// RKB-assembled relativistic kinetic (Orbital_RKB_IBS_Imp::MakeKinetic), which already folds in the
// c_light factors, so NO extra 1/2 or c is applied here.  CAUTION: those c_light/2 factors in the RKB
// chain (atom Orbital_RKBL_IBS::MakeKinetic, Imp/Orbital_DHF_IBS.C) appear to cancel but are
// UNVERIFIED -- to be sorted in a dedicated relativistic-kinetic cleanup.  See BasisSet/Orbital_1E_IBS.C.
rsmat_t DiracKinetic::CalculateMatrix(const obs_t* bs,const Spin&) const
{
    return bs->Kinetic();
}

void DiracKinetic::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    te.Kinetic=cd->DM_Contract(this);
}

std::ostream& DiracKinetic::Write(std::ostream& os) const
{
    os << "    Dirac kinetic energy c*sigma*p" << std::endl;
    return os;
}

} //namespace

