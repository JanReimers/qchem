#include <gtkmm.h>
#include <iostream>

#include "Imp/Cluster/Atom.H"
#include "Imp/BasisSet/Atom/l/Slater_BS.H"
#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/ml/Gaussian_BS.H"
#include "Imp/Hamiltonian/Hamiltonians.H"

class AtomFrame : public Gtk::Frame
{
public:
  AtomFrame();
  AtomFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~AtomFrame();

  Atom* create() const;

private:
  const Glib::RefPtr<Gtk::Builder> itsBuilder;
  Gtk::SpinButton* itsZ_spin;
  Gtk::SpinButton* itsCharge_spin;
};

AtomFrame::AtomFrame() : Glib::ObjectBase("atom_frame")
{
}

AtomFrame::AtomFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder)
: Glib::ObjectBase("atom_frame"),
  Gtk::Frame(cobject),
  itsBuilder(refBuilder),
  itsZ_spin(refBuilder->get_widget<Gtk::SpinButton>("z_spin")),
  itsCharge_spin(refBuilder->get_widget<Gtk::SpinButton>("charge_spin"))
  {
    
  }
  
AtomFrame::~AtomFrame()
{ 
}
Atom* AtomFrame::create() const
{
    int Z=itsZ_spin->get_value_as_int();
    int charge=itsCharge_spin->get_value_as_int();
    return new Atom(Z,charge);
}

class BasisSetFrame : public Gtk::Frame
{
public:
  BasisSetFrame();
  BasisSetFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~BasisSetFrame();

  BasisSet* create() const;

private:
  enum bstypes {SlaterYl,SlaterYlm,GaussianYl,GaussianYlm};
  static const std::map<Glib::ustring,bstypes> bstype_map;
  static bstypes find(Glib::ustring);

  const Glib::RefPtr<Gtk::Builder> itsBuilder;
  Gtk::DropDown* itsType;
  Glib::RefPtr<Gtk::StringList> itsTypes; 
  Gtk::Entry* itsEmin;
  Gtk::Entry* itsEmax;
  Gtk::SpinButton* itsN;
};

const std::map<Glib::ustring,BasisSetFrame::bstypes> BasisSetFrame::bstype_map=
  {
    {"Slater Yl",SlaterYl},
    {"Slater Ylm",SlaterYlm},
    {"Gaussian Yl",GaussianYl},
    {"Gaussian Ylm",GaussianYlm},
  };

BasisSetFrame::BasisSetFrame() : Glib::ObjectBase("basisset_frame")
{}
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

BasisSet* BasisSetFrame::create() const
{
  guint it=itsType->get_selected();
  Glib::ustring bs_stype=itsTypes->get_string(it);  
  bstypes bs_type=find(bs_stype);
  double emin=Glib::Ascii::strtod(itsEmin->get_text());
  double emax=Glib::Ascii::strtod(itsEmax->get_text());
  guint N=itsN->get_value_as_int();
  size_t LMax=0;
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

class HamiltonianFrame : public Gtk::Frame
{
public:
  typedef std::shared_ptr<const Cluster> cl_t;
  HamiltonianFrame();
  HamiltonianFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~HamiltonianFrame();

  Hamiltonian* create(const cl_t& cl,const MeshParams& m, const BasisSet* bs) const;

private:
  enum htypes {H1E,HF,DFT,D1E,DHF};
  static const std::map<Glib::ustring,htypes> htype_map;
  static htypes find(Glib::ustring);

  const Glib::RefPtr<Gtk::Builder> itsBuilder;
  Gtk::DropDown* itsType;
  Glib::RefPtr<Gtk::StringList> itsTypes; 
  Gtk::CheckButton* itsPolarized;
};

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
  
Hamiltonian* HamiltonianFrame::create(const cl_t& cl,const MeshParams& m, const BasisSet* bs) const
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
      h= polarized ? (Hamiltonian*)new Ham_DFT_P(cl,0.7,m,bs) : (Hamiltonian*)new Ham_DFT_U(cl,0.7,m,bs);
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


class Controller
{
public:
  Controller(const Glib::RefPtr<Gtk::Builder>& refBuilder)
  : itsStartButton(refBuilder->get_widget<Gtk::Button>("start"))
  , itsStepButton(refBuilder->get_widget<Gtk::Button>("step"))
  , itsStopButton(refBuilder->get_widget<Gtk::Button>("pause"))
  , itsAtom(Gtk::Builder::get_widget_derived<AtomFrame>(refBuilder, "atom_frame"))
  , itsBasisSet(Gtk::Builder::get_widget_derived<BasisSetFrame>(refBuilder, "basisset_frame"))
  , itsHamiltonian(Gtk::Builder::get_widget_derived<HamiltonianFrame>(refBuilder, "ham_frame"))
  {};

private:
  Gtk::Button* itsStartButton;
  Gtk::Button* itsStepButton;
  Gtk::Button* itsStopButton;
  AtomFrame* itsAtom;
  BasisSetFrame* itsBasisSet;
  HamiltonianFrame* itsHamiltonian;
};







Glib::RefPtr<Gtk::Application> app;
   
void on_app_activate()
{
    auto refBuilder = Gtk::Builder::create();
    try
    {
      refBuilder->add_from_file("qchem-gtk.ui");
    }
    catch(const Glib::FileError& ex)
    {
      std::cerr << "FileError: " << ex.what() << std::endl;
      return;
    }
    catch(const Glib::MarkupError& ex)
    {
      std::cerr << "MarkupError: " << ex.what() << std::endl;
      return;
    }
    catch(const Gtk::BuilderError& ex)
    {
      std::cerr << "BuilderError: " << ex.what() << std::endl;
      return;
    }
  
    Controller c(refBuilder);
    auto main = refBuilder->get_widget<Gtk::Window>("main");
    app->add_window(*main);
    main->set_visible(true);
}

int main(int argc, char** argv)
{
    app = Gtk::Application::create("org.gtkmm.example");

    app->signal_activate().connect([] () { on_app_activate(); });

  return app->run(argc, argv);
}