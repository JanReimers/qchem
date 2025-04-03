// File: MasterUnPolarizedWF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/MasterUnPolarizedWF.H"
#include "Imp/SCFIterator/SCFIteratorUnPol.H"
#include "Imp/WaveFunction/IrrepWaveFunction.H"
#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/Orbitals/TOrbitals.H"
#include <BasisSet.H>
#include <Irrep_BS.H>
// #include <Symmetry.H>
#include <cassert>
#include <Spin.H>
#include <Orbital_QNs.H>
#include <iostream>

MasterUnPolarizedWF::MasterUnPolarizedWF()
    : itsEC(0)
{};

MasterUnPolarizedWF::MasterUnPolarizedWF(const BasisSet* bs,const ElectronConfiguration* ec)
    : itsBS(bs)
    , itsEC(ec)
{
    assert(itsBS);
    assert(itsEC);
    assert(itsBS->GetNumFunctions()>0);
    for (auto b:itsBS->Iterate<TOrbital_IBS<double> >())
    {
        uiwf_t wfp(new IrrepWaveFunction(b,Spin(Spin::None)));
        itsIWFs.push_back(std::move(wfp));
    }
};

MasterUnPolarizedWF::~MasterUnPolarizedWF()
{
    
}
//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void MasterUnPolarizedWF::DoSCFIteration(Hamiltonian& ham)
{
    for (auto& w:itsIWFs) w->DoSCFIteration(ham);
}

Exact_CD* MasterUnPolarizedWF::GetChargeDensity(Spin s) const
{
    assert(s==Spin::None);
    // return itsGroup->GetChargeDensity(s);
    
    Composite_Exact_CD* cd = new Composite_Exact_CD();
    for (auto& w:itsIWFs) cd->Insert(w->GetChargeDensity(s));
    return cd;
}

Orbitals* MasterUnPolarizedWF::GetOrbitals(const Irrep_QNs& qns) const
{
    // No UT coverage
    assert(qns.ms==Spin::None);
    auto i=itsQN_WFs.find(qns);
    assert(i!=itsQN_WFs.end());
    return i->second->GetOrbitals(qns);
}

SCFIterator* MasterUnPolarizedWF::MakeIterator(Hamiltonian* H, Exact_CD* cd, double nElectrons)
{
    return new SCFIteratorUnPol(this, H, cd,nElectrons);
}

const EnergyLevels& MasterUnPolarizedWF::FillOrbitals(const ElectronConfiguration*)
{
    itsELevels.clear();
    for (auto& w:itsIWFs) 
        itsELevels.merge(w->FillOrbitals(itsEC),0.0001);
    return itsELevels;
}

void MasterUnPolarizedWF::DisplayEigen() const
{
    std::cout << "Alpha+Beta spin :" << std::endl;
    itsELevels.Report(std::cout);
}


std::ostream& MasterUnPolarizedWF::Write(std::ostream& os) const
{
    return os;
}

std::istream& MasterUnPolarizedWF::Read (std::istream& is)
{
    return is;
}

