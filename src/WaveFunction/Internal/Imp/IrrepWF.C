// File: IrrepWF.C  Wave function for an unpolarized atom.
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
#include "blaze/Math.h" 
import qchem.SCFAccelerator;
module qchem.WaveFunction.Internal.IrrepWF;
import qchem.Orbitals.Factory;
import qchem.Blaze;

namespace qchem::WaveFunction
{

using std::cout;
using std::endl;

IrrepWF::IrrepWF(const obs_t* bs, LASolver<double>* lasb,const Irrep& qns,SCFIrrepAccelerator* acc)
    : itsBasisSet   (bs)
    , itsLASolver   (lasb)
    , itsOrbitals   (qchem::Orbitals::Factory(bs,qns.ms))
    , itsIrrep      (qns)
    , itsAccelerator(acc)
    , itsDPrime     (zero<double>(bs->GetNumFunctions()))
{
    assert(itsOrbitals);
    assert(itsAccelerator);
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
    // The accelerator returns the next orbitals: DIIS extrapolates F' and diagonalizes;
    // a direct minimizer (GDM) rotates the current orbitals along the Grassmann manifold.
    auto [U,Up,e]=itsAccelerator->NextOrbitals();
    itsOrbitals->UpdateOrbitals(U,Up,e);
}

// Direct-minimization: ask the accelerator to compute its step (false in the seed step).
bool IrrepWF::ComputeStep()
{
    assert(itsOrbitals);
    return itsAccelerator->ComputeStep();
}
// Direct-minimization: move the orbitals to geodesic fraction t (commit=false is a trial).
void IrrepWF::MoveOrbitals(double t, bool commit)
{
    assert(itsOrbitals);
    auto [U,Up,e]=itsAccelerator->OrbitalsAt(t,commit);
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
    for (auto o:itsOrbitals->Iterate<qchem::Orbitals::Orbital>())
        itsELevels.insert(qchem::Orbitals::EnergyLevel(o));
    
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

rvec_t IrrepWF::Get_BS_Diagonal() const
{
    assert(itsLASolver);
    return itsLASolver->Get_BS_Diagonal();
}

void  IrrepWF::DisplayEigen() const
{
    itsELevels.Report(std::cout);
   
}

} //namespace
