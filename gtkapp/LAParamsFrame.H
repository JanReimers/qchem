// File: LAParamsFrame.H  GTK frame to show and manage linear algebra settings.
#ifndef _LAParamsFrame_H_
#define _LAParamsFrame_H_

#include <gtkmm.h>
#include <types.H>
#include <LAParams.H>

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

#endif //_LAParamsFrame_H_