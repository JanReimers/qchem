// File: HamiltonianFrame.cpp  GTK frame to show and manage Hamiltonian settings.

#include "HamiltonianFrame.H"
#include "Hamiltonians.H"
#include <iostream>


HamiltonianFrame::HamiltonianFrame() {};
HamiltonianFrame::HamiltonianFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder)
  : Glib::ObjectBase("ham_frame")
  , Gtk::Frame(cobject)
  , itsEnumDD(Gtk::Builder::get_widget_derived<enumDropDown<htypes>>(refBuilder, "ham_dropdown"))
  , itsPolarized(refBuilder->get_widget<Gtk::CheckButton>("ham_polarized"))
{
  itsEnumDD->init({H1E,HF,DFT,D1E,DHF},{"1-Electron (1E)","Hatree-Fock (HF)","Density-Functional (DFT)","Dirac 1E","Dirac HF"});
  
}

HamiltonianFrame::~HamiltonianFrame() {};


 
Hamiltonian* HamiltonianFrame::create(const cl_t& cl,const MeshParams* m, const BasisSet* bs) const
{
  h_type=itsEnumDD->GetType();
  is_polarized=itsPolarized->get_active();
  Hamiltonian* h=0;
  switch (h_type)
  {
    case H1E : 
      h=new Ham_1E(cl);
      break;
    case HF : 
      h= is_polarized ? (Hamiltonian*)new Ham_HF_P(cl) : (Hamiltonian*)new Ham_HF_U(cl);
      break;
    case DFT : 
      assert(m);
      h= is_polarized ? (Hamiltonian*)new Ham_DFT_P(cl,0.7,*m,bs) : (Hamiltonian*)new Ham_DFT_U(cl,0.7,*m,bs);
      break;
    case D1E : 
      h= new Ham_DHF_1E(cl);
      break;
    case DHF : 
      h= new Ham_DHF(cl);
      break;
  } 
  return h;
}

#include "Imp/WaveFunction/UnPolarized_WF.H"
#include "Imp/WaveFunction/Polarized_WF.H"

WaveFunction* HamiltonianFrame::create(BasisSet* bs, ElectronConfiguration* ec ) const
{
    assert(bs);
    assert(ec);
    if (itsPolarized->get_active())
        return new Polarized_WF(bs,ec);
    else
        return new UnPolarized_WF(bs,ec);
}

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