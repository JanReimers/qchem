// File: IrrepWF.C  Wave function for an unpolarized atom.
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
import qchem.SCFAccelerator;
module qchem.WaveFunction.Internal.IrrepWF;
import qchem.Orbitals.Factory;

using std::cout;
using std::endl;

IrrepWF::IrrepWF(const TOrbital_IBS<double>* bs, LASolver<double>* las, const Irrep_QNs& qns,SCFIrrepAccelerator* acc)
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

IrrepWF::~IrrepWF()
{
    delete itsOrbitals;
    delete itsLASolver;
    delete itsAccelerator;
}

void IrrepWF::CalculateH(Hamiltonian& ham,const DM_CD* cd)
{
    assert(itsOrbitals);
    itsF=ham.GetMatrix(itsBasisSet,itsIrrep.ms,cd); //Hamiltonian or Fock matrix in the non-orthogonal basis.
    itsAccelerator->UseFD(itsF,itsDPrime); //Feed non-ortho F into the accelerator along with density matrix (in the orthogonal basis).
}

//----------------------------------------------------------------------------
//
//  This function will create unoccupied orbtials.  
//
void IrrepWF::DoSCFIteration()
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
const EnergyLevels& IrrepWF::FillOrbitals(const ElectronConfiguration* ec)
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


DM_CD* IrrepWF::GetChargeDensity() const
{
    assert(itsOrbitals);
    return itsOrbitals->GetChargeDensity();
}

const Orbitals* IrrepWF::GetOrbitals() const
{
    assert(itsOrbitals);
    return itsOrbitals;
}
Orbitals* IrrepWF::GetOrbitals() 
{
    assert(itsOrbitals);
    return itsOrbitals;
}

Vector<double> IrrepWF::Get_BS_Diagonal() const
{
    assert(itsLASolver);
    return itsLASolver->Get_BS_Diagonal();
}

void  IrrepWF::DisplayEigen() const
{
    itsELevels.Report(std::cout);
   
}

