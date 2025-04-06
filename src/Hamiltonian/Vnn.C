// File: Vnn.C  Nuclear-Nuclear potential.



#include "Imp/Cluster/Atom.H"
#include "Imp/Hamiltonian/Vnn.H"
#include <TotalEnergy.H>
#include <Irrep_BS.H>
#include "oml/smatrix.h"
#include "oml/vector.h"
#include "oml/vector3d.h"
#include <cassert>

Vnn::Vnn()
    : Static_HT_Imp()
    , theCluster()
{};

Vnn::Vnn(cl_t& cl)
    : Static_HT_Imp()
    , theCluster(cl)
{};

Static_HT::SMat Vnn::CalculateHamiltonianMatrix(const ibs_t* bs,const Spin&) const
{
    int n=bs->GetNumFunctions();
    SMat ret(n,n);
    Fill(ret,0.0); //No contribution to Hamiltonian matrix here.
    return ret;
}

void Vnn::GetEnergy(TotalEnergy& te,const Exact_CD* cd) const
{
    double vnn=0.0;
    for(auto atom1:*theCluster)
        for(auto atom2:*theCluster)
        {
            RVec3 r1=atom1->itsR, r2=atom2->itsR;
            if (r1!=r2) vnn += 0.5 * atom1->itsZ * atom2->itsZ / norm(r1-r2);
        }

    te.Enn=vnn;
}

std::ostream& Vnn::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        size_t Na=theCluster->GetNumAtoms();
        os << "    Nuclear-Nuclear potential ZiZj/|Ri-Rj| with " << Na*(Na-1) << " nucleus pairs." << std::endl;
    }
    else
        os << *theCluster;
    return os;
}

std::istream& Vnn::Read (std::istream& is)
{
    return is;
}

