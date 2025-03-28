// File: HamiltonianImplementation.C  General matrix implementation of a Hamiltonian operator.



#include "Imp/Hamiltonian/Hamiltonian.H"
#include "Imp/Hamiltonian/Kinetic.H"
#include "Imp/Hamiltonian/Ven.H"
#include "Imp/Hamiltonian/Vnn.H"
#include <TotalEnergy.H>
#include <ChargeDensity.H>
#include <Irrep_BS.H>
#include "Imp/Containers/ptr_vector_io.h"
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>

typedef optr_vector1<HamiltonianTerm*>::iterator        ITER;
typedef optr_vector1<HamiltonianTerm*>::const_iterator CITER;

HamiltonianImp::HamiltonianImp()
    : itsExactCD(0)
{};

void HamiltonianImp::Add(HamiltonianTerm* p)
{
    itsHamiltonianTerms.push_back(p);
}

void HamiltonianImp::InsertStandardTerms(cl_t & cl)
{
    Add(new Kinetic);
    Add(new Vnn(cl));
    Add(new Ven(cl));
}

void HamiltonianImp::UseChargeDensity(const Exact_CD* cd)
{
    itsExactCD =cd;
    assert(itsExactCD);
    for (auto t:itsHamiltonianTerms) t->UseChargeDensity(itsExactCD);
}

//
//  If any of the HamiltonianTerm terms is polarized, then Hamiltonian is considered polarized.
//
bool HamiltonianImp::IsPolarized() const
{
    //No UT coverage
    bool ret=false;
    for (auto t:itsHamiltonianTerms) ret = ret || t->IsPolarized();
    return ret;
}

Hamiltonian::SMat HamiltonianImp::BuildHamiltonian(const TOrbital_IBS<double>* bs,const Spin& S) const
{
    int n=bs->GetNumFunctions();
    SMat H(n,n);
    Fill(H,0.0);
    for (auto t:itsHamiltonianTerms) H+=t->BuildHamiltonian(bs,S);
    return H;
}


TotalEnergy HamiltonianImp::GetTotalEnergy() const
{
    assert(itsExactCD);
    TotalEnergy e;
    for (auto t:itsHamiltonianTerms)  t->GetEnergy(e);
    return e;
}


std::ostream& HamiltonianImp::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "Hamiltonian with " << itsHamiltonianTerms.size() << " terms:" << std::endl;
    os << itsHamiltonianTerms;
    return os;
}

std::istream& HamiltonianImp::Read(std::istream& is)
{
    is >> itsHamiltonianTerms;
    return is;
}


