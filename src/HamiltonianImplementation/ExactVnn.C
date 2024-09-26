// File: ExactVnn.C  Nuclear potential.



#include "HamiltonianImplementation/ExactVnn.H"
#include "Hamiltonian/TotalEnergy.H"
#include "Cluster/ClusterBrowser.H"
#include "BasisSet/BasisSet.H"
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
    for(ClusterBrowser b1(*theCluster); b1; b1++)
        for(ClusterBrowser b2(*theCluster); b2; b2++)
        {
            RVec3 r1=(*b1).itsR, r2=(*b2).itsR;
            if (r1!=r2) vnn += 0.5 * (*b1).itsZ * (*b2).itsZ / !(r1-r2);
        }

    te.Enn=vnn;
}

std::ostream& ExactVnn::Write(std::ostream& os) const
{
    return os << *theCluster;
}

std::istream& ExactVnn::Read (std::istream& is)
{
    theCluster.reset(Cluster::Factory(is));
    return is >> *theCluster;
}

