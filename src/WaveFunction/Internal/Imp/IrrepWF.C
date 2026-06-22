// File: IrrepWF.C  Wave function for an unpolarized atom.
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
module qchem.WaveFunction.Internal.IrrepWF;
import qchem.SCFAccelerator;
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
    , itsDPrime     (blazem::zero<double>(bs->GetNumFunctions()))
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
    return FillOrbitals((double)ec->GetN(itsIrrep)); // electrons for this Irrep={spin,symmetry}
}

const EnergyLevels& IrrepWF::FillOrbitals(double ne)
{
    // Empty first: the aufbau occupation can shift between iterations, and TakeElectrons only
    // overwrites the orbitals it fills (leaving stale occupation on the rest).
    for (auto o:itsOrbitals->Iterate()) o->Empty();
    std::tie(ne,itsDPrime)=itsOrbitals->TakeElectrons(ne); // occupy lowest-first, build density
    assert(ne==0.0); //enough orbitals to take all electrons; if not the basis set is too small.

    // List of energy levels.  Degenerate levels should get merged.
    itsELevels.clear();
    for (auto o:itsOrbitals->Iterate())
        itsELevels.insert(qchem::Orbitals::EnergyLevel(o));

    return itsELevels;
}


// MOM: score each current orbital by the norm of its projection onto the reference occupied
// subspace.  Coefficients are in the orthonormal basis (metric = I), so the overlap of orbital j
// with reference orbital i is just the dot product; the score is sqrt(sum_i <ref_i|j>^2) in [0,1].
// Storage order matches Iterate<Orbital> (same underlying vector), so scores align with the aufbau.
rvec_t IrrepWF::MOMScores() const
{
    assert(itsOrbitals);
    if (itsRefOccCPrime.columns()==0) return rvec_t();          // no reference captured yet
    std::vector<double> sc;
    for (auto o:itsOrbitals->Iterate<qchem::Orbitals::TOrbital<double>>())
    {
        rvec_t proj=blazem::trans(itsRefOccCPrime)*o->GetCoeffPrime();   // (nref) overlaps
        sc.push_back(blazem::norm(proj));
    }
    rvec_t scores(sc.size());
    for (size_t i=0;i<sc.size();++i) scores[i]=sc[i];
    return scores;
}

// Snapshot the currently-occupied orbitals' C' columns as the reference for the next iteration.
// An unoccupied irrep clears its reference (no continuation to track).
void IrrepWF::CaptureMOMReference()
{
    assert(itsOrbitals);
    std::vector<vec_t<double>> cols;
    for (auto o:itsOrbitals->Iterate<qchem::Orbitals::TOrbital<double>>())
        if (o->IsOccupied()) cols.push_back(o->GetCoeffPrime());
    if (cols.empty()) { itsRefOccCPrime.clear(); return; }
    itsRefOccCPrime.resize(cols.front().size(),cols.size());
    for (size_t j=0;j<cols.size();++j) blazem::column(itsRefOccCPrime,j)=cols[j];
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
