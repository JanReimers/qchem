// File: HamiltonianImplementation.C  General matrix implementation of a Hamiltonian operator.



#include "TotalEnergy.H"
#include "HamiltonianImplementation/HamiltonianImplementation.H"
#include "ChargeDensity.H"
#include "BasisSet.H"
#include "Imp/Containers/ptr_vector_io.h"
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>

typedef optr_vector1<HamiltonianTerm*>::iterator        ITER;
typedef optr_vector1<HamiltonianTerm*>::const_iterator CITER;

HamiltonianImplementation::HamiltonianImplementation()
    : itsExactCD(0)
{};

void HamiltonianImplementation::Add(HamiltonianTerm* p)
{
    itsHamiltonianTerms.push_back(p);
}

void HamiltonianImplementation::UseChargeDensity(const ChargeDensity* theExactCD)
{
    itsExactCD =theExactCD;
    assert(itsExactCD);
    for (auto t:itsHamiltonianTerms) t->UseChargeDensity(itsExactCD);
//    for (ITER i(itsHamiltonianTerms.begin()); i!=itsHamiltonianTerms.end(); i++)
//        i->UseChargeDensity(itsExactCD);
}

//
//  If any of the HamiltonianTerm terms is polarized, then Hamiltonian is considered polarized.
//
bool HamiltonianImplementation::IsPolarized() const
{
    //No UT coverage
    bool ret=false;
    for (auto t:itsHamiltonianTerms) ret = ret || t->IsPolarized();
    return ret;
}

Hamiltonian::SMat HamiltonianImplementation::BuildHamiltonian(const IrrepBasisSet* bs,const Spin& S) const
{
    int n=bs->GetNumFunctions();
    SMat H(n,n);
    Fill(H,0.0);
    for (auto t:itsHamiltonianTerms) H+=t->BuildHamiltonian(bs,S);
    return H;
}


TotalEnergy HamiltonianImplementation::GetTotalEnergy() const
{
    assert(itsExactCD);
    TotalEnergy e;
    for (auto t:itsHamiltonianTerms)  t->GetEnergy(e);
    return e;
}


std::ostream& HamiltonianImplementation::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "Hamiltonian with " << itsHamiltonianTerms.size() << " terms:" << std::endl;
    os << itsHamiltonianTerms;
    return os;
}

std::istream& HamiltonianImplementation::Read(std::istream& is)
{
    is >> itsHamiltonianTerms;
    return is;
}


