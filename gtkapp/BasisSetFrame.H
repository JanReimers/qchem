// File: BasisSetFrame.H  GTK frame to show and manage basis set settings.
#ifndef _BasisSetFrame_H_
#define _BasisSetFrame_H_

#include <gtkmm.h>
#include <types.H>

class BasisSetFrame : public Gtk::Frame
{
public:
  BasisSetFrame();
  BasisSetFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~BasisSetFrame();

  BasisSet* create(size_t LMax) const;

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


#endif //_BasisSetFrame_H_