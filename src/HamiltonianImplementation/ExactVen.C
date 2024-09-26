// File: ExactVen.C  Nuclear potential.



#include "HamiltonianImplementation/ExactVen.H"
#include "BasisSet/BasisSet.H"
#include "Hamiltonian/TotalEnergy.H"
#include "ChargeDensity/ChargeDensity.H"
#include "oml/smatrix.h"
#include "oml/vector.h"
#include "oml/vector3d.h"

ExactVen::ExactVen()
    : HamiltonianTermImplementation()
    , theCluster()
{};

ExactVen::ExactVen(const rc_ptr<Cluster>& cl)
    : HamiltonianTermImplementation()
    , theCluster(cl)
{
    assert(cl->GetNumAtoms()>0);
};

//double ExactVen::operator()(const RVec3& r) const
//{
//    return 1.0/!r;
//}
//
//ExactVen::Vec3 ExactVen::Gradient(const RVec3& r) const
//{
//    return -r/!r;
//}

HamiltonianTerm::SMat ExactVen::CalculateHamiltonianMatrix(const BasisSet* bs,const Spin&) const
{
    return bs->GetNuclear(&*theCluster);
}

void ExactVen::GetEnergy(TotalEnergy& te) const
{
    te.Een=CalculateEnergy();
}

std::ostream& ExactVen::Write(std::ostream& os) const
{
    return os << *theCluster;
}

std::istream& ExactVen::Read (std::istream& is)
{
    Cluster* cl=Cluster::Factory(is);
    is >> *cl;
    theCluster.reset(cl);
    return is;
}

