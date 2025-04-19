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
 
    itsNotebook->set_margin(10);
    itsNotebook->set_expand();

    // auto plotw=new PlotWindow();
    // auto plotw=new Gtk::Label("asdfasd");
    // itsNotebook->append_page(*plotw,"plot1");

};

#include "Imp/Cluster/Molecule.H"
#include "Imp/Cluster/Atom.H"
#include <MeshParams.H>
#include <SCFIterator.H>
#include <IterationParams.H>
#include <BasisSet.H>

void ControllerWindow::new_model()
{
  std::cout << "New model!!" << std::endl;
  LAParams lap=itsLAParams->create();

  HamiltonianFrame::cl_t a(itsAtom->create());
  BasisSet* bs=itsBasisSet->create();
  bs->Set(lap);
  MeshParams m=itsMeshFrame->create();
  Hamiltonian* h=itsHamiltonian->create(a,&m,bs);
  Atom_EC ec(a->GetNumElectrons());
  WaveFunction* wf=itsHamiltonian->create(bs,&ec);
  SCFIterationParams scfip=itsIterationParams->create();
  itsSCFIterator=new SCFIterator(wf,h);
  itsSCFIterator->Iterate(scfip);
  auto plotw=itsHamiltonian->create_pw(bs,wf);
  itsNotebook->append_page(*plotw,"Orbitals");
  
}


