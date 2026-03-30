// File: Vnn.C  Nuclear-Nuclear potential.
module;
#include <iostream>
#include <cassert>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.IrrepBasisSet;
import qchem.Atom;
import qchem.Blaze;

Vnn::Vnn(const cl_t& cl)
    : Static_HT_Imp()
    , theCluster(cl)
{};

rsmat_t Vnn::CalculateMatrix(const ibs_t* bs,const Spin&) const
{
    int n=bs->GetNumFunctions();
    rsmat_t ret=zero<double>(n);
    return ret;
}

void Vnn::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    double vnn=0.0;
    for(const auto& atom1:*theCluster)
        for(const auto& atom2:*theCluster)
        {
            Vector3D<double> r1=atom1->itsR, r2=atom2->itsR;
            if (r1!=r2) vnn += 0.5 * atom1->itsZ * atom2->itsZ / norm(r1-r2);
        }

    te.Enn=vnn;
}

std::ostream& Vnn::Write(std::ostream& os) const
{
    size_t Na=theCluster->GetNumAtoms();
    os << "    Nuclear-Nuclear potential ZiZj/|Ri-Rj| with " << Na*(Na-1) << " nucleus pairs." << std::endl;
    return os;
}

