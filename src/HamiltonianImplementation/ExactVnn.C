// File: ExactVnn.C  Nuclear potential.



#include "HamiltonianImplementation/ExactVnn.H"
#include "TotalEnergy.H"
#include "BasisSet.H"
#include "oml/smatrix.h"
#include "oml/vector.h"
#include "oml/vector3d.h"
#include <cassert>

ExactVnn::ExactVnn()
    : HamiltonianTermImplementation()
    , theCluster()
{};

ExactVnn::ExactVnn(const rc_ptr<Cluster>& cl)
    : HamiltonianTermImplementation()
    , theCluster(cl)
{};

HamiltonianTerm::SMat ExactVnn::CalculateHamiltonianMatrix(const BasisSet* bs,const Spin&) const
{
    int n=bs->GetNumFunctions();
    SMat ret(n,n);
    Fill(ret,0.0); //No contribution to Hamiltonian matrix here.
    return ret;
}

void ExactVnn::GetEnergy(TotalEnergy& te) const
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

std::ostream& ExactVnn::Write(std::ostream& os) const
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

std::istream& ExactVnn::Read (std::istream& is)
{
    theCluster.reset(Cluster::Factory(is));
    return is >> *theCluster;
}

