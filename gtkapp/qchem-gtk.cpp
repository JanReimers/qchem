#include <gtkmm.h>
#include <iostream>

#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
#include "Imp/BasisSet/Atom/l/Slater_BS.H"
#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/ml/Gaussian_BS.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Containers/stl_io.h"
#include "Imp/Mesh/LogRadialMesh.H"
#include "Imp/Mesh/MHLRadialMesh.H"
#include "Imp/Mesh/GaussAngularMesh.H"
#include "Imp/Mesh/GaussLegendreAngularMesh.H"
#include "Imp/Mesh/EulerMaclarenAngularMesh.H"
#include "Imp/Cluster/AtomMesh.H"


#include <LAParams.H>
#include <IterationParams.H>
#include <MeshParams.H>

class AtomFrame : public Gtk::Frame
{
public:
  AtomFrame();
  AtomFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~AtomFrame();

  Molecule* create() const;

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
Molecule* AtomFrame::create() const
{
    int Z=itsZ_spin->get_value_as_int();
    int charge=itsCharge_spin->get_value_as_int();
    Molecule* m=new Molecule();
    m->Insert(new Atom(Z,charge));
    return m;
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

  Hamiltonian* create(const cl_t& cl,const MeshParams* m, const BasisSet* bs) const;

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

class LAParamsFrame : public Gtk::Frame
{
public:
  LAParamsFrame() {};
  LAParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~LAParamsFrame() {};

  LAParams create() const;
private:
  
  static const std::map<Glib::ustring,qchem::Pkg> pkgtype_map;
  static const std::map<Glib::ustring,qchem::Ortho> orthotype_map;
  static qchem::Pkg   find_pkg(Glib::ustring);
  static qchem::Ortho find_ortho(Glib::ustring);

  Glib::RefPtr<Gtk::StringList> itsPkgTypes; 
  Glib::RefPtr<Gtk::StringList> itsOrthoTypes; 

  Gtk::DropDown* itsPkgType;
  Gtk::DropDown* itsOrthoType;
  Gtk::Entry*    itsTruncationTolerance;
  Gtk::Entry*    itsAbsTol;
};

const std::map<Glib::ustring,qchem::Pkg> LAParamsFrame::pkgtype_map=
{
  {"LAPack",qchem::Lapack},
  {"OML"   ,qchem::OML},
};

const std::map<Glib::ustring,qchem::Ortho> LAParamsFrame::orthotype_map=
{
  {"Cholsky",qchem::Cholsky},
  {"SVD"    ,qchem::SVD},
  {"Eigen"  ,qchem::Eigen}, 
};

LAParamsFrame::LAParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("laparams_frame")
  , Gtk::Frame(cobject)
  , itsPkgType(refBuilder->get_widget<Gtk::DropDown>("la_pkg"))
  , itsOrthoType(refBuilder->get_widget<Gtk::DropDown>("la_ortho"))
  , itsTruncationTolerance(refBuilder->get_widget<Gtk::Entry>("la_trunc"))
  , itsAbsTol(refBuilder->get_widget<Gtk::Entry>("la_abstol"))
{
  std::vector<Glib::ustring> strings;
  for (const auto& [key, _] : pkgtype_map) strings.push_back(key);
  itsPkgTypes = Gtk::StringList::create(strings);
  itsPkgType->set_model(itsPkgTypes);
  itsPkgType->set_selected(0);

  strings.clear();
  for (const auto& [key, _] : orthotype_map) strings.push_back(key);
  itsOrthoTypes = Gtk::StringList::create(strings);
  itsOrthoType->set_model(itsOrthoTypes);
  itsOrthoType->set_selected(0);
}

qchem::Pkg LAParamsFrame::find_pkg(Glib::ustring s)
{
    auto i=pkgtype_map.find(s);
    if (i==pkgtype_map.end())
    {
      std::cerr << "LAParamsFrame::find Unknown package type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

qchem::Ortho LAParamsFrame::find_ortho(Glib::ustring s)
{
    auto i=orthotype_map.find(s);
    if (i==orthotype_map.end())
    {
      std::cerr << "LAParamsFrame::find Unknown ortho type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

LAParams LAParamsFrame::create() const
{
  guint it=itsPkgType->get_selected();
  Glib::ustring pkg_stype=itsPkgTypes->get_string(it);  
  qchem::Pkg  pkg_type=find_pkg(pkg_stype);
  it=itsOrthoType->get_selected();
  Glib::ustring ortho_stype=itsOrthoTypes->get_string(it);  
  qchem::Ortho  ortho_type=find_ortho(ortho_stype);

  double trunctol=Glib::Ascii::strtod(itsTruncationTolerance->get_text());
  double abstol=Glib::Ascii::strtod(itsAbsTol->get_text());
  return {pkg_type,ortho_type,trunctol,abstol};
}

class SCFIterationParamsFrame : public Gtk::Frame
{
public:
SCFIterationParamsFrame() {};
SCFIterationParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~SCFIterationParamsFrame() {};

  SCFIterationParams create() const;
private:
    
  Gtk::SpinButton*  itsNIter;
  Gtk::Entry*       itsMinDeltaRo;
  Gtk::Entry*       itsRoRelax;
  Gtk::CheckButton* itsVerbose;
};

SCFIterationParamsFrame::SCFIterationParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("iteration_frame")
  , Gtk::Frame(cobject)
  , itsNIter(refBuilder->get_widget<Gtk::SpinButton>("SCF_max_iter"))
  , itsMinDeltaRo(refBuilder->get_widget<Gtk::Entry>("SCF_delta_ro"))
  , itsRoRelax(refBuilder->get_widget<Gtk::Entry>("SCF_ro_relax"))
  , itsVerbose(refBuilder->get_widget<Gtk::CheckButton>("SCF_verbose"))
{
  
}

SCFIterationParams SCFIterationParamsFrame::create() const
{
  guint N_iter=itsNIter->get_value_as_int();
  double min_delta_ro=Glib::Ascii::strtod(itsMinDeltaRo->get_text());
  double ro_relax=Glib::Ascii::strtod(itsRoRelax->get_text());
  bool verbose=itsVerbose->get_active();
  return {N_iter,min_delta_ro,ro_relax,0.0,verbose};
}

class MeshFrame : public Gtk::Frame
{
public:
  MeshFrame() {};
  MeshFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~MeshFrame() {};

  Mesh* create() const;
private:
  
  static const std::map<Glib::ustring,qchem::RadialType> radial_type_map;
  static const std::map<Glib::ustring,qchem::AngleType> angular_type_map;
  static qchem::RadialType find_radial(Glib::ustring);
  static qchem::AngleType  find_angular(Glib::ustring);

  Glib::RefPtr<Gtk::StringList> itsRadialTypes; 
  Glib::RefPtr<Gtk::StringList> itsAngularTypes; 

  Gtk::DropDown* itsRadialType;
  Gtk::DropDown* itsAngularType;
  Gtk::SpinButton* itsNRadial;
  Gtk::SpinButton* itsNAngular;
};

const std::map<Glib::ustring,qchem::RadialType> MeshFrame::radial_type_map=
{
  {"MHL",qchem::MHL},
  {"Log",qchem::Log},
};
const std::map<Glib::ustring,qchem::AngleType> MeshFrame::angular_type_map=
{
  {"Gauss"         ,qchem::Gauss},
  {"Gauss Legendre",qchem::GaussLegendre},
  {"Euler Mclaren" ,qchem::EulerMclaren},
};

MeshFrame::MeshFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("mesh_frame")
  , Gtk::Frame(cobject)
  , itsRadialType (refBuilder->get_widget<Gtk::DropDown>("mesh_radial"))
  , itsAngularType(refBuilder->get_widget<Gtk::DropDown>("mesh_angular"))
  , itsNRadial    (refBuilder->get_widget<Gtk::SpinButton>("mesh_nr"))
  , itsNAngular   (refBuilder->get_widget<Gtk::SpinButton>("mesh_na"))
{
  std::vector<Glib::ustring> strings;
  for (const auto& [key, _] : radial_type_map) strings.push_back(key);
  itsRadialTypes = Gtk::StringList::create(strings);
  itsRadialType->set_model(itsRadialTypes);
  itsRadialType->set_selected(0);

  strings.clear();
  for (const auto& [key, _] : angular_type_map) strings.push_back(key);
  itsAngularTypes = Gtk::StringList::create(strings);
  itsAngularType->set_model(itsAngularTypes);
  itsAngularType->set_selected(0);
}

qchem::RadialType MeshFrame::find_radial(Glib::ustring s)
{
    auto i=radial_type_map.find(s);
    if (i==radial_type_map.end())
    {
      std::cerr << "MeshFrame::find Unknown radial mesh type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

qchem::AngleType MeshFrame::find_angular(Glib::ustring s)
{
    auto i=angular_type_map.find(s);
    if (i==angular_type_map.end())
    {
      std::cerr << "MeshFrame::find Unknown angular mesh type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

Mesh* MeshFrame::create() const
{
  guint it=itsRadialType->get_selected();
  Glib::ustring r_stype=itsRadialTypes->get_string(it);  
  qchem::RadialType r_type=find_radial(r_stype);
  it=itsAngularType->get_selected();
  Glib::ustring a_stype=itsAngularTypes->get_string(it);  
  qchem::AngleType a_type=find_angular(a_stype);

  guint Nr=itsNRadial->get_value_as_int();
  guint Na=itsNAngular->get_value_as_int();
  int m=2,L=Na;
  double alpha=2.0,start=0.01,stop=40.0;

  RadialMesh* mr=0;
  switch (r_type)
  {
    case qchem::MHL : 
      mr=new MHLRadialMesh(Nr,m,alpha);
      break;
    case qchem::Log : 
      mr=new LogRadialMesh(start,stop,Nr);
      break;
  } 

  Mesh* ma=0;
  switch (a_type)
  {
    case qchem::Gauss :
      ma= new GaussAngularMesh(Na);
      break;
    case qchem::GaussLegendre :
      ma= new GaussLegendreAngularMesh(L,m);
      break;
    case qchem::EulerMclaren :
      ma= new EulerMaclarenAngularMesh(L,m);
  }
  return new AtomMesh(*mr,*ma,RVec3(0,0,0));
}


class ControllerWindow : public Gtk::Window
{
public:
  ControllerWindow(BaseObjectType* cobject,const Glib::RefPtr<Gtk::Builder>& refBuilder)
  : Glib::ObjectBase("main")
  , Gtk::Window(cobject)
  , itsStartButton(refBuilder->get_widget<Gtk::Button>("start"))
  , itsStepButton(refBuilder->get_widget<Gtk::Button>("step"))
  , itsPauseButton(refBuilder->get_widget<Gtk::Button>("pause"))
  , itsAtom(Gtk::Builder::get_widget_derived<AtomFrame>(refBuilder, "atom_frame"))
  , itsBasisSet(Gtk::Builder::get_widget_derived<BasisSetFrame>(refBuilder, "basisset_frame"))
  , itsHamiltonian(Gtk::Builder::get_widget_derived<HamiltonianFrame>(refBuilder, "ham_frame"))
  , itsLAParams(Gtk::Builder::get_widget_derived<LAParamsFrame>(refBuilder, "laparams_frame"))
  , itsIterationParams(Gtk::Builder::get_widget_derived<SCFIterationParamsFrame>(refBuilder, "Iteration_frame"))
  , itsMeshFrame(Gtk::Builder::get_widget_derived<MeshFrame>(refBuilder, "mesh_frame"))
  {
    itsStartButton->signal_clicked().connect(sigc::mem_fun(*this,&ControllerWindow::new_model));
    itsStepButton->signal_clicked().connect(sigc::mem_fun(*this,&ControllerWindow::new_model));
  };

  void new_model();

private:
  Gtk::Button* itsStartButton;
  Gtk::Button* itsStepButton;
  Gtk::Button* itsPauseButton;
  AtomFrame*        itsAtom;
  BasisSetFrame*    itsBasisSet;
  HamiltonianFrame* itsHamiltonian;
  LAParamsFrame*    itsLAParams;
  SCFIterationParamsFrame* itsIterationParams;
  MeshFrame* itsMeshFrame;
};


void ControllerWindow::new_model()
{
  std::cout << "New model!!" << std::endl;
  HamiltonianFrame::cl_t a(itsAtom->create());
  BasisSet* bs=itsBasisSet->create();
  Hamiltonian* h=itsHamiltonian->create(a,0,bs);
}




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
  
    // Controller c(refBuilder);
    // auto main = refBuilder->get_widget<Gtk::Window>("main");
    auto main= Gtk::Builder::get_widget_derived<ControllerWindow>(refBuilder, "main");
    app->add_window(*main);
    main->set_visible(true);
}

int main(int argc, char** argv)
{
    app = Gtk::Application::create("org.gtkmm.example");

    app->signal_activate().connect([] () { on_app_activate(); });

  return app->run(argc, argv);
}