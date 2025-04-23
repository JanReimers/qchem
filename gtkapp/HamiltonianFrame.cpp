// File: HamiltonianFrame.cpp  GTK frame to show and manage Hamiltonian settings.

#include "HamiltonianFrame.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include <iostream>


HamiltonianFrame::HamiltonianFrame() {};
HamiltonianFrame::HamiltonianFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder)
  : Glib::ObjectBase("ham_frame")
  , Gtk::Frame(cobject)
  , itsBuilder(refBuilder)
  , itsType(refBuilder->get_widget<Gtk::DropDown>("ham_dropdown"))
  , itsPolarized(refBuilder->get_widget<Gtk::CheckButton>("ham_polarized"))
{
  std::vector<Glib::ustring> strings;
  for (const auto& [key, _] : htype_map) strings.push_back(key);

  itsTypes = Gtk::StringList::create(strings); //This auto sorts alphabetically.
  htype_invmap.clear();
  for (guint i=0;i<itsTypes->get_n_items();i++)
  {
    htypes ht=find(itsTypes->get_string(i));
    htype_invmap[ht]=i;
  }
  itsType->set_model(itsTypes);
  itsType->set_selected(0);
}

HamiltonianFrame::~HamiltonianFrame() {};

void HamiltonianFrame::init()
{
  // std::cout << "h_type=" << h_type << " " << htype_invmap[h_type] << std::endl;
  itsType->set_selected(htype_invmap[h_type]); //implicit conversion of enum to int .. one off problem?
  itsPolarized->set_active(is_polarized);
}

const std::map<Glib::ustring,HamiltonianFrame::htypes> HamiltonianFrame::htype_map=
{ 
{"1-Electron (1E)",H1E},
{"Hatree-Fock (HF)",HF},
{"Density-Functional (DFT)",DFT},
{"Dirac 1E",D1E},
{"Dirac HF",DHF},
};

std::map<HamiltonianFrame::htypes,guint> HamiltonianFrame::htype_invmap;

HamiltonianFrame::htypes HamiltonianFrame::find(Glib::ustring s)
{
    auto i=htype_map.find(s);
    if (i==htype_map.end())
    {
        std::cerr << "HamiltonianFrame::find Unknown Hamiltonian type '" << s << "'" << std::endl;
        exit(-1);
    }
    return i->second;
}
  
HamiltonianFrame::htypes HamiltonianFrame::GetType() const
{
  guint it=itsType->get_selected();
  Glib::ustring h_stype=itsTypes->get_string(it);  
  return find(h_stype);
}

Hamiltonian* HamiltonianFrame::create(const cl_t& cl,const MeshParams* m, const BasisSet* bs) const
{
  h_type=GetType();
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