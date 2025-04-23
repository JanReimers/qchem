// File: ControllerWindow.cpp  GTK main window for the atom GUI app.

#include "ControllerWindow.H"
#include "AtomFrame.H"
#include "BasisSetFrame.H"
#include "HamiltonianFrame.H"
#include "LAParamsFrame.H"
#include "SCFFrame.H"
#include "MeshFrame.H"
#include "Imp/BasisSet/Atom/EC.H"
#include "PlotWindow.H"
#include <fstream>
#include <cereal/archives/json.hpp>

std::string ControllerWindow::defaultFilename=getenv("HOME")+std::string("/.qchem");

ControllerWindow::ControllerWindow(BaseObjectType* cobject,const Glib::RefPtr<Gtk::Builder>& refBuilder)
: Glib::ObjectBase("main")
, Gtk::Window(cobject)
, itsNotebook(refBuilder->get_widget<Gtk::Notebook>("notebook"))
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
  this->signal_close_request().connect(sigc::mem_fun(*this,&ControllerWindow::close_request),false);
 
    itsNotebook->set_margin(10);
    itsNotebook->set_expand();

  {
    
    std::ifstream fs(defaultFilename);
    if (fs)
    {
      std::cout << "Reading archive" << std::endl;
      cereal::JSONInputArchive iarchive(fs); 
      iarchive(*this);
    }
  }
};


ControllerWindow::~ControllerWindow()
{
  // std::cout << "~ControllerWindow" << std::endl;
}

bool ControllerWindow::close_request()
{
  // std::cout << "Dumping archive: '" << defaultFilename << "'" << std::endl;
  std::ofstream fs(defaultFilename);
  cereal::JSONOutputArchive oarchive(fs); 
  oarchive(*this);
  return false;
}

#include "Imp/Cluster/Molecule.H"
#include "Imp/Cluster/Atom.H"
#include <MeshParams.H>
#include <SCFIterator.H>
#include <IterationParams.H>
#include <BasisSet.H>
#include <WaveFunction.H>
#include <ChargeDensity.H>


void ControllerWindow::new_model()
{
  std::cout << "New model!!" << std::endl;
  for (int i=0;i<itsNotebook->get_n_pages();i++)
    itsNotebook->remove_page(i);
  LAParams lap=itsLAParams->create();

  HamiltonianFrame::cl_t a(itsAtom->create());
  Atom_EC ec(a->GetNumElectrons());
  BasisSet* bs=itsBasisSet->create(ec.GetLMax());
  bs->Set(lap);
  MeshParams m=itsMeshFrame->create();
  Hamiltonian* h=itsHamiltonian->create(a,&m,bs);
  WaveFunction* wf=itsHamiltonian->create(bs,&ec);
  SCFIterationParams scfip=itsIterationParams->create();
  itsSCFIterator=new SCFIterator(wf,h);
  itsSCFIterator->Iterate(scfip);
  auto orbital_plotw=itsHamiltonian->create_orbital_pw(bs,wf);
  orbital_plotw->AddLabels("#fir (a.u.)","#fi#gf(r)","Radial Orbital Profiles");
  itsNotebook->append_page(*orbital_plotw,"φ(r)");

  typedef ScalarFunction<double> sf_t;

  sf_t* cd=wf->GetChargeDensity();
  sf_t* sd=wf->GetSpinDensity();
  itsNotebook->append_page(*new SF_PW(cd,"#fir (a.u.)","#fir#u2#d ρ(r) (a.u.)","Radial Charge Density"),"ρ(r)");
  if (sd) itsNotebook->append_page(*new SF_PW(sd,"#fir (a.u.)","#fir#u2#d(ρ#d↑#u-ρ#d↓#u) (a.u.)","Radial Spin Density"),"ρ↑-ρ↓");
  
}


