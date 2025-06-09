// File: Irrep_WF.C  Wave function for an unpolarized atom.

#include "Imp/WaveFunction/Irrep_WF.H"
#include "Imp/SCFAccelerator.H"
#include "Imp/Orbitals/TOrbitals.H"

#include <ElectronConfiguration.H>
#include <Irrep_BS.H>
#include <LASolver.H>
#include <Hamiltonian.H>
#include "oml/vector.h"

using std::cout;
using std::endl;

Irrep_WF::Irrep_WF(const TOrbital_IBS<double>* bs, const Spin& ms,SCFIrrepAccelerator* acc)
    : itsBasisSet(bs)
    , itsLASolver(bs->CreateSolver())
    , itsOrbitals(new  TOrbitalsImp<double>(bs,ms))
    , itsIrrep     (ms,bs->GetSymmetry())
    , itsAccelerator(acc)
    , itsDPrime(bs->GetNumFunctions())
{
    assert(itsOrbitals);
    assert(itsAccelerator);
    itsAccelerator->Init(itsLASolver,itsIrrep);
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
    itsF=ham.GetMatrix(itsBasisSet,itsIrrep.ms,cd);
    itsAccelerator->UseFD(itsF,itsDPrime);
}

//----------------------------------------------------------------------------
//
//  This function will create unoccupied orbtials.  
//
void Irrep_WF::DoSCFIteration(Hamiltonian& ham,const DM_CD* cd)
{
    assert(itsOrbitals);
    // itsF=ham.GetMatrix(itsBasisSet,itsIrrep.ms,cd);
    //
    //  Feed F,D into the SCF accelerator
    //
    SMatrix<double> FPrime=itsAccelerator->Project(); //Calcularte F'=Vd*F*V and conditionally project F'
    // cout << "FPrime=" << FPrime << endl;
    auto [U,Up,e]=itsLASolver->SolveOrtho(FPrime);
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

