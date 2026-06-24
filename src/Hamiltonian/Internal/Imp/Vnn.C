// File: Vnn.C  Nuclear-Nuclear potential.
module;
#include <iostream>
#include <cassert>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.Structure;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

Vnn::Vnn(const st_t& st)
    : Static_HT_Imp()
    , theStructure(st)
{};

rsmat_t Vnn::CalculateMatrix(const obs_t* bs,const Spin&) const
{
    int n=bs->GetNumFunctions();
    rsmat_t ret=blazem::zero<double>(n);
    return ret;
}

void Vnn::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    double vnn=0.0;
    for(const auto& atom1:*theStructure)
        for(const auto& atom2:*theStructure)
        {
            Vector3D<double> r1=atom1->itsR, r2=atom2->itsR;
            if (r1!=r2) vnn += 0.5 * atom1->itsZ * atom2->itsZ / norm(r1-r2);
        }

    te.Enn=vnn;
}

std::ostream& Vnn::Write(std::ostream& os) const
{
    size_t Na=theStructure->GetNumAtoms();
    os << "    Nuclear-Nuclear potential ZiZj/|Ri-Rj| with " << Na*(Na-1) << " nucleus pairs." << std::endl;
    return os;
}

} //namespace
