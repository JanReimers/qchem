// File: MeshFrame.H  GTK frame to show and manage DFT integration mesh settings.
#ifndef _MeshFrame_H_
#define _MeshFrame_H_

#include <gtkmm.h>
#include <types.H>
#include <MeshParams.H>

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

#endif //_MeshFrame_H_