// File: HamiltonianImplementation.C  General matrix implementation of a Hamiltonian operator.


#include "Common/stl_io.h"
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>

#include "Hamiltonian.H"
#include "Kinetic.H"
#include "Ven.H"
#include "Vnn.H"
#include <Hamiltonian/TotalEnergy.H>
#include <ChargeDensity/ChargeDensity.H>
#include <BasisSet/Irrep_BS.H>

HamiltonianImp::HamiltonianImp() : itsIsPolarized(false)
{};

void HamiltonianImp::Add(Static_HT* p)
{
    itsSHTs.push_back(std::unique_ptr<Static_HT>(p));
    itsIsPolarized = itsIsPolarized || p->IsPolarized();
}
void HamiltonianImp::Add(Dynamic_HT* p)
{
    itsDHTs.push_back(std::unique_ptr<Dynamic_HT>(p));
    itsIsPolarized = itsIsPolarized || p->IsPolarized();
}

void HamiltonianImp::InsertStandardTerms(const cl_t & cl)
{
    Add(new Kinetic);
    Add(new Vnn(cl));
    Add(new Ven(cl));
}

Hamiltonian::SMat HamiltonianImp::GetMatrix(const ibs_t* bs,const Spin& S,const DM_CD* cd)
{
    int n=bs->GetNumFunctions();
    SMat H(n,n);
    Fill(H,0.0);
    for (auto& t:itsSHTs) H+=t->GetMatrix(bs,S);
    // Leave these terms out if we don't have guess for the charge density.
    if (cd)
        for (auto& t:itsDHTs) H+=t->GetMatrix(bs,S,cd);
    return H;
}


EnergyBreakdown HamiltonianImp::GetTotalEnergy( const DM_CD* cd ) const
{
    assert(cd);
    EnergyBreakdown e;
    for (auto& t:itsSHTs)  t->GetEnergy(e,cd);
    for (auto& t:itsDHTs)  t->GetEnergy(e,cd);
    return e;
}


std::ostream& HamiltonianImp::Write(std::ostream& os) const
{
    if (itsIsPolarized) os << "Polarized ";
    os << "Hamiltonian with " << itsSHTs.size() << " static terms:" << std::endl;
    os << itsSHTs;
    if (itsIsPolarized) os << "Polarized ";
    os << "Hamiltonian with " << itsDHTs.size() << " dynamic terms:" << std::endl;
    os << itsDHTs;
    return os;
}
