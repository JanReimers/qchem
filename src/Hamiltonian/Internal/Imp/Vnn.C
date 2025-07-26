// File: Vnn.C  Nuclear-Nuclear potential.
module;
#include <iostream>
#include <cassert>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.Irrep_BS;
import qchem.Atom;

Vnn::Vnn()
    : Static_HT_Imp()
    , theCluster()
{};

Vnn::Vnn(const cl_t& cl)
    : Static_HT_Imp()
    , theCluster(cl)
{};

 SMatrix<double>  Vnn::CalculateMatrix(const ibs_t* bs,const Spin&) const
{
    int n=bs->GetNumFunctions();
    SMatrix<double> ret(n,n);
    Fill(ret,0.0); //No contribution to Hamiltonian matrix here.
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

std::istream& Vnn::Read (std::istream& is)
{
    return is;
}

