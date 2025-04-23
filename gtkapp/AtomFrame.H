// File: AtomFrame.H GTK frame to show and manage atom settings.
#ifndef _AtomFrame_H_
#define _AtomFrame_H_

#include <gtkmm.h>
#include <types.H>

class AtomFrame : public Gtk::Frame
{
public:
  AtomFrame();
  AtomFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~AtomFrame();

  Molecule* create() const;
  template<class Archive>  void save(Archive & archive) const
  {
    Z=itsZ_spin->get_value_as_int();
    charge=itsCharge_spin->get_value_as_int();
    archive( Z,charge); // serialize things by passing them to the archive
  }
  template<class Archive>  void load(Archive & archive)
  {
    archive( Z,charge); 
    itsZ_spin->set_value(Z);
    itsCharge_spin->set_value(charge);
  }
private:
  const Glib::RefPtr<Gtk::Builder> itsBuilder;
  Gtk::SpinButton* itsZ_spin;
  Gtk::SpinButton* itsCharge_spin;

  mutable int Z,charge;
};

#endif //_AtomFrame_H_
