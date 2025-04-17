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

  itsTypes = Gtk::StringList::create(strings);
  itsType->set_model(itsTypes);
  itsType->set_selected(0);
}

HamiltonianFrame::~HamiltonianFrame() {};

const std::map<Glib::ustring,HamiltonianFrame::htypes> HamiltonianFrame::htype_map=
{
{"1-Electron (1E)",H1E},
{"Hatree-Fock (HF)",HF},
{"Density-Functional (DFT)",DFT},
{"Dirac 1E",D1E},
{"Dirac HF",DHF},
};

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
  
Hamiltonian* HamiltonianFrame::create(const cl_t& cl,const MeshParams* m, const BasisSet* bs) const
{
  guint it=itsType->get_selected();
  Glib::ustring h_stype=itsTypes->get_string(it);  
  htypes h_type=find(h_stype);
  bool polarized=itsPolarized->get_active();
  Hamiltonian* h=0;
  switch (h_type)
  {
    case H1E : 
      h=new Ham_1E(cl);
      break;
    case HF : 
      h= polarized ? (Hamiltonian*)new Ham_HF_P(cl) : (Hamiltonian*)new Ham_HF_U(cl);
      break;
    case DFT : 
      assert(m);
      h= polarized ? (Hamiltonian*)new Ham_DFT_P(cl,0.7,*m,bs) : (Hamiltonian*)new Ham_DFT_U(cl,0.7,*m,bs);
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

