// File: BasisSetFrame.cpp  GTK frame to show and manage basis set settings.

#include "BasisSetFrame.H"
#include "Imp/BasisSet/Atom/l/Slater_BS.H"
#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/ml/Gaussian_BS.H"

const std::map<Glib::ustring,BasisSetFrame::bstypes> BasisSetFrame::bstype_map=
  {
    {"Slater Yl",SlaterYl},
    {"Slater Ylm",SlaterYlm},
    {"Gaussian Yl",GaussianYl},
    {"Gaussian Ylm",GaussianYlm},
  };

BasisSetFrame::BasisSetFrame() : Glib::ObjectBase("basisset_frame") {}

BasisSetFrame::BasisSetFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("basisset_frame")
  , Gtk::Frame(cobject)
  , itsBuilder(refBuilder)
  , itsType(refBuilder->get_widget<Gtk::DropDown>("basisset_dropdown"))
  , itsEmin(refBuilder->get_widget<Gtk::Entry>("basisset_emin"))
  , itsEmax(refBuilder->get_widget<Gtk::Entry>("basisset_emax"))
  , itsN(refBuilder->get_widget<Gtk::SpinButton>("basisset_N"))
{
  std::vector<Glib::ustring> strings;
  for (const auto& [key, _] : bstype_map) strings.push_back(key);
  
  itsTypes = Gtk::StringList::create(strings);
  itsType->set_model(itsTypes);
  itsType->set_selected(0);
}
  
BasisSetFrame::~BasisSetFrame() {};

BasisSetFrame::bstypes BasisSetFrame::find(Glib::ustring s)
{
    auto i=bstype_map.find(s);
    if (i==bstype_map.end())
    {
      std::cerr << "BasisSetFrame::find Unknown basis set type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

BasisSet* BasisSetFrame::create(size_t LMax) const
{
  guint it=itsType->get_selected();
  Glib::ustring bs_stype=itsTypes->get_string(it);  
  bstypes bs_type=find(bs_stype);
  double emin=Glib::Ascii::strtod(itsEmin->get_text());
  double emax=Glib::Ascii::strtod(itsEmax->get_text());
  guint N=itsN->get_value_as_int();
  BasisSet* bs=0;
  switch (bs_type)
  {
    case SlaterYl : 
      bs=new Atoml::Slater::BasisSet(N,emin,emax,LMax);
      break;
    case SlaterYlm : 
      bs=new Atom_ml::Slater::BasisSet(N,emin,emax,LMax);
      break;
    case GaussianYl : 
      bs= new Atoml::Gaussian::BasisSet(N,emin,emax,LMax);
      break;
    case GaussianYlm : 
      bs= new Atom_ml::Gaussian::BasisSet(N,emin,emax,LMax);
      break;
  } 
  return bs;
}
