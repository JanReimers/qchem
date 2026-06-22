// File: AtomFrame.cpp GTK frame to show and manage atom settings.

#include "AtomFrame.H"
import qchem.Structure.Atom;

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
Structure* AtomFrame::create() const
{
    Z=itsZ_spin->get_value_as_int();
    charge=itsCharge_spin->get_value_as_int();
    return new Atom(Z,charge);
}

