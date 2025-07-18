// File: Irrep_WF.C  Wave function for an unpolarized atom.

#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
#include "Irrep_WF.H"
#include <SCFAccelerator/SCFAccelerator.H>

#include <Symmetry/ElectronConfiguration.H>
#include <BasisSet/Irrep_BS.H>
#include <LASolver/LASolver.H>
#include <Hamiltonian/Hamiltonian.H>
#include <Orbitals/Factory.H>
#include <Orbitals/Orbitals.H>

using std::cout;
using std::endl;

Irrep_WF::Irrep_WF(const TOrbital_IBS<double>* bs, LASolver<double>* las, const Irrep_QNs& qns,SCFIrrepAccelerator* acc)
    : itsBasisSet   (bs)
    , itsLASolver   (las)
    , itsOrbitals   (OrbitalsF::Factory(bs,qns.ms))
    , itsIrrep      (qns)
    , itsAccelerator(acc)
    , itsDPrime     (bs->GetNumFunctions())
{
    assert(itsOrbitals);
    assert(itsAccelerator);
    Fill(itsDPrime,0.0);
};

Irrep_WF::~Irrep_WF()
{
    delete itsOrbitals;
    delete itsLASolver;
    delete itsAccelerator;
}

void Irrep_WF::CalculateH(Hamiltonian& ham,const DM_CD* cd)
{
    assert(itsOrbitals);
    itsF=ham.GetMatrix(itsBasisSet,itsIrrep.ms,cd); //Hamiltonian or Fock matrix in the non-orthogonal basis.
    itsAccelerator->UseFD(itsF,itsDPrime); //Feed non-ortho F into the accelerator along with density matrix (in the orthogonal basis).
}

//----------------------------------------------------------------------------
//
//  This function will create unoccupied orbtials.  
//
void Irrep_WF::DoSCFIteration()
{
    assert(itsOrbitals);
    //project F' using pre calculated coefficients. And then diagonalize it.
    auto [U,Up,e]=itsLASolver->SolveOrtho(itsAccelerator->Project());
    itsOrbitals->UpdateOrbitals(U,Up,e);
}
//
//  Now populate the orbitals with electrons.  The ElectronConfiguration knows how many electrons
//  are in each Irrep.
//
const EnergyLevels& Irrep_WF::FillOrbitals(const ElectronConfiguration* ec)
{
    
    double ne=ec->GetN(itsIrrep); // Step one: How many electron for this Irrep={spin,symmetry} ?
    std::tie(ne,itsDPrime)=itsOrbitals->TakeElectrons(ne); // Step two: Dump electrons into the orbitals and then calculate a density matrix.
    assert(ne==0.0); //There must be enough orbitals to take all electrons for the Irrep.  If not the basis set is too small.
    
    // Step three: Make a list of energy levels.  Degenerate levels should get merged.
    itsELevels.clear();
    for (auto o:itsOrbitals->Iterate<Orbital>())
        itsELevels.insert(EnergyLevel(o));
    
    return itsELevels;
}


DM_CD* Irrep_WF::GetChargeDensity() const
{
    assert(itsOrbitals);
    return itsOrbitals->GetChargeDensity();
}

const Orbitals* Irrep_WF::GetOrbitals() const
{
    assert(itsOrbitals);
    return itsOrbitals;
}
Orbitals* Irrep_WF::GetOrbitals() 
{
    assert(itsOrbitals);
    return itsOrbitals;
}

Vector<double> Irrep_WF::Get_BS_Diagonal() const
{
    assert(itsLASolver);
    return itsLASolver->Get_BS_Diagonal();
}

void  Irrep_WF::DisplayEigen() const
{
    itsELevels.Report(std::cout);
   
}

