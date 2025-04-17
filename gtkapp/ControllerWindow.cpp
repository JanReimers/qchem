// File: ControllerWindow.cpp  GTK main window for the atom GUI app.

#include "ControllerWindow.H"
#include "AtomFrame.H"
#include "BasisSetFrame.H"
#include "HamiltonianFrame.H"
#include "LAParamsFrame.H"
#include "SCFFrame.H"
#include "MeshFrame.H"

ControllerWindow::ControllerWindow(BaseObjectType* cobject,const Glib::RefPtr<Gtk::Builder>& refBuilder)
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

#include "Imp/Cluster/Molecule.H"
#include "Imp/Cluster/Atom.H"
void ControllerWindow::new_model()
{
  std::cout << "New model!!" << std::endl;
  HamiltonianFrame::cl_t a(itsAtom->create());
  BasisSet* bs=itsBasisSet->create();
  //Hamiltonian* h=
  itsHamiltonian->create(a,0,bs);
}


