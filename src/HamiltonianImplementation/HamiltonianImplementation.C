// File: HamiltonianImplementation.C  General matrix implementation of a Hamiltonian operator.



#include "Hamiltonian/TotalEnergy.H"
#include "HamiltonianImplementation/HamiltonianImplementation.H"
#include "ChargeDensity/ChargeDensity.H"
#include "BasisSet/BasisSet.H"
#include "Misc/ptrvector_io.h"
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>

typedef optr_vector<HamiltonianTerm*>::iterator        ITER;
typedef optr_vector<HamiltonianTerm*>::const_iterator CITER;

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
    for (ITER i(itsHamiltonianTerms.begin()); i!=itsHamiltonianTerms.end(); i++)
        i->UseChargeDensity(itsExactCD);
}

//
//  If any of the HamiltonianTerm terms is polarized, then Hamiltonian is considered polarized.
//
bool HamiltonianImplementation::IsPolarized() const
{
    bool ret=false;
    for (CITER b(itsHamiltonianTerms.begin()); b!=itsHamiltonianTerms.end(); b++)
        ret = ret || b->IsPolarized();
    return ret;
}

Hamiltonian::SMat HamiltonianImplementation::BuildHamiltonian(const BasisSet* bs,const Spin& S) const
{
    int n=bs->GetNumFunctions();
    SMat H(n,n);
    Fill(H,0.0);
    for(CITER b(itsHamiltonianTerms.begin()); b!=itsHamiltonianTerms.end(); b++)
        H+=b->BuildHamiltonian(bs,S);
    return H;
}


TotalEnergy HamiltonianImplementation::GetTotalEnergy() const
{
    assert(itsExactCD);
    TotalEnergy e;
    for(CITER p(itsHamiltonianTerms.begin()); p!=itsHamiltonianTerms.end(); p++) p->GetEnergy(e);
    return e;
}


std::ostream& HamiltonianImplementation::Write(std::ostream& os) const
{
    os << itsHamiltonianTerms;
    return os;
}

std::istream& HamiltonianImplementation::Read(std::istream& is)
{
    is >> itsHamiltonianTerms;
    return is;
}


