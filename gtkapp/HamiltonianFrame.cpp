// File: HamiltonianFrame.cpp  GTK frame to show and manage Hamiltonian settings.

#include "HamiltonianFrame.H"
#include <iostream>

import qchem.Hamiltonians;
using namespace qchem::Hamiltonian;

HamiltonianFrame::HamiltonianFrame() {};
HamiltonianFrame::HamiltonianFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder)
  : Glib::ObjectBase("ham_frame")
  , Gtk::Frame(cobject)
  , itsEnumDD(Gtk::Builder::get_widget_derived<enumDropDown<htypes>>(refBuilder, "ham_dropdown"))
  , itsPolarized(refBuilder->get_widget<Gtk::CheckButton>("ham_polarized"))
{
  itsEnumDD->init({E1,HF,DE1,DHF},{"1-Electron (1E)","Hatree-Fock (HF)","Dirac 1E","Dirac HF"});
  
}

HamiltonianFrame::~HamiltonianFrame() {};


 
Hamiltonian* HamiltonianFrame::create(const cl_t& cl,const MeshParams* m, const BasisSet* bs) const
{
  h_type=itsEnumDD->GetType();
  is_polarized=itsPolarized->get_active();
  Pol pol= itsPolarized->get_active() ? Pol::Polarized : Pol::UnPolarized;
  return Factory(h_type,pol,cl);
}

// #include "Imp/WaveFunction/UnPolarizedWF.H"
// #include "Imp/WaveFunction/PolarizedWF.H"

// WaveFunction* HamiltonianFrame::create(BasisSet* bs, ElectronConfiguration* ec ) const
// {
//     assert(bs);
//     assert(ec);
//     if (itsPolarized->get_active())
//         return new PolarizedWF(bs,ec);
//     else
//         return new UnPolarizedWF(bs,ec);
// }

#include "PlotWindow.H"
PlotWindow* HamiltonianFrame::create_orbital_pw(BasisSet* bs,WaveFunction* wf) const
{
  assert(bs);
  assert(wf);
  if (itsPolarized->get_active())
        return new Polarized_Orbital_PW(bs,wf);
    else
        return new UnPolarized_Orbital_PW(bs,wf);
}

PlotWindow* HamiltonianFrame::create_el_pw(WaveFunction* wf) const
{
  assert(wf);
  if (itsPolarized->get_active())
        return new Polarized_EnergyLevel_PW(wf);
    else
        return new EnergyLevel_PW(wf);
}

PlotWindow*   HamiltonianFrame::create_diag_pw(BasisSet* bs,WaveFunction* wf, qchem::Ortho ortho) const
{
  Spin s = itsPolarized->get_active() ? Spin::Up : Spin::None;
  return new Diagonal_PW(bs,wf,ortho,s);
}