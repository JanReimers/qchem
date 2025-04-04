// File: MasterPolarizedWF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/MasterPolarizedWF.H"
#include "Imp/SCFIterator/SCFIteratorPol.H"
#include "Imp/ChargeDensity/PolarizedCD.H"
#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/WaveFunction/IrrepWaveFunction.H"
#include "Imp/Orbitals/TOrbitals.H"
#include <BasisSet.H>
#include <Symmetry.H>
#include <Irrep_BS.H>
#include <cassert>
#include <Spin.H>
#include <Orbital_QNs.H>
#include <iostream>



MasterPolarizedWF::MasterPolarizedWF(const BasisSet* bs,const ElectronConfiguration* ec)
    : itsBS(bs) //Basis set
    , itsEC(ec) //Electron cofiguration
{
    assert(itsEC);
    assert(itsBS->GetNumFunctions()>0);
    for (auto b:itsBS->Iterate<TOrbital_IBS<double> >())
    {
        uiwf_t wfup(new IrrepWaveFunction(b,Spin(Spin::Up)));
        uiwf_t wfdn(new IrrepWaveFunction(b,Spin(Spin::Down)));
        itsQN_WFs[wfup->GetQNs()]=wfup.get();
        itsQN_WFs[wfdn->GetQNs()]=wfdn.get();
        // Do tte move last.
        itsSpinUpIWFs.push_back(std::move(wfup));
        itsSpinDnIWFs.push_back(std::move(wfdn));
    }
};

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void MasterPolarizedWF::DoSCFIteration(Hamiltonian& ham)
{
    for (auto& w:itsSpinUpIWFs) w->DoSCFIteration(ham);
    for (auto& w:itsSpinDnIWFs) w->DoSCFIteration(ham);
}

Exact_CD* MasterPolarizedWF::GetChargeDensity(Spin s) const
{
    Composite_Exact_CD* up = new Composite_Exact_CD();
    Composite_Exact_CD* dn = new Composite_Exact_CD();
    for (auto& w:itsSpinUpIWFs) up->Insert(w->GetChargeDensity(s));
    for (auto& w:itsSpinDnIWFs) dn->Insert(w->GetChargeDensity(s));
    return new Polarized_CDImp(up,dn);
}

Orbitals* MasterPolarizedWF::GetOrbitals(const Irrep_QNs& qns) const
{
    auto i=itsQN_WFs.find(qns);
    assert(i!=itsQN_WFs.end());
    return i->second->GetOrbitals(qns);
}


const EnergyLevels& MasterPolarizedWF::FillOrbitals(const ElectronConfiguration*)
{
    itsUpELevels.clear();
    itsDnELevels.clear();

    for (auto& w:itsSpinUpIWFs) 
         itsUpELevels.merge(w->FillOrbitals(itsEC),0.0001);
    for (auto& w:itsSpinDnIWFs) 
         itsDnELevels.merge(w->FillOrbitals(itsEC),0.0001);
    return itsUpELevels;
}


SCFIterator* MasterPolarizedWF::MakeIterator(Hamiltonian* H, Exact_CD* cd, double nElectrons)
{
    return new SCFIteratorPol(this,H,cd,nElectrons,0.0);
}

void MasterPolarizedWF::DisplayEigen() const
{
    std::cout << "Spin:         up                 down            avg" << std::endl;
    auto iup=itsUpELevels.begin();
    auto idn=itsDnELevels.begin();
    while (iup!=itsUpELevels.end() && idn!=itsDnELevels.end())
    {
        bool valid_up = iup!=itsUpELevels.end();
        bool valid_dn = idn!=itsDnELevels.end();
        bool print_up = valid_up && iup->first<=0.0;
        bool print_dn = valid_dn && idn->first<=0.0;
        bool qn_match = valid_up && valid_dn && *iup->second.qn==*idn->second.qn;
        
        if (print_up)
            iup->second.Report(std::cout);
        else
            std::cout << "                                     ";
        
        if (print_dn && qn_match)
            idn->second.Report(std::cout);
        
        
            
        std::cout << std::endl;
        if (valid_up) iup++;
        if (valid_dn && qn_match) idn++; //If no match, let iup catch up.
        if (iup->first>0.0 && idn->first>0.0) break;
    }
   
}

std::ostream& MasterPolarizedWF::Write(std::ostream& os) const
{
    return os;
}

std::istream& MasterPolarizedWF::Read (std::istream& is)
{
    return is;
}


