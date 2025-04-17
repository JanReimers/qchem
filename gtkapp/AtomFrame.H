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

private:
  const Glib::RefPtr<Gtk::Builder> itsBuilder;
  Gtk::SpinButton* itsZ_spin;
  Gtk::SpinButton* itsCharge_spin;
};

#endif //_AtomFrame_H_
