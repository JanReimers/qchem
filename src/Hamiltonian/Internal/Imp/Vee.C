// File: Vee.C  Electron-Electron Coulomb potential
module;
#include <iostream>
#include <vector>
#include <map>
#include <string>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Types;
import qchem.ChargeDensity;
import qchem.Energy;

namespace qchem::Hamiltonian
{

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro(r2)/r12 .
//              /
//  Where ro is the charge density.
//

// Whole-system Coulomb (doc/ERI4Rework.md §5.4): the density scatters itself across canonical irrep pairs
// (AccumulateDirectAll -> ScatterBoth), so only the canonical ERI4 blocks are ever built/cached -- halving
// Coulomb ERI RAM+build vs the old per-irrep pass, which fetched both J(i,j) and J(j,i).  Result is bit-close
// (~1e-13) to the old Fock, not bit-identical: the old J(j,i) was an independent build, here it is J(i,j)^T.
void Vee::AccumulateAll(std::vector<rsmat_t>& X,const std::vector<const ohfbs_t*>& hf,const DM_CD* dm) const
{
    dm->AccumulateDirectAll(X,hf);
}

void Vee::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    // E_ee = 1/2 Tr(D.J), taken from THIS term's own whole-system Coulomb blocks -- no per-irrep GetMatrix
    // round-trip through DM_Contract (which is what kept Vee tied to the 3-arg GetMatrix / tDynamic_CC).
    ContractAll(cd);
    te.Eee=0.5*cd->DM_ContractBlocks(itsJKs);
    te.EeeFit    = 0.0;
    te.EeeFitFit = 0.0;
}

std::ostream& Vee::Write(std::ostream& os) const
{
    os << "    Coulomb potential ro(r_2)/r_12" << std::endl;
    return os;
}

} //namespace
