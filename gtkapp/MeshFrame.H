// File: MeshFrame.H  GTK frame to show and manage DFT integration mesh settings.
#ifndef _MeshFrame_H_
#define _MeshFrame_H_

#include "enumDropDown.H"
#include <gtkmm.h>
#include <types.H>
#include <Mesh/MeshParams.H>

class MeshFrame : public Gtk::Frame
{
public:
  MeshFrame() {};
  MeshFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~MeshFrame() {};

  MeshParams create() const;

  template<class Archive>  void save(Archive & archive) const
  {
    rtype=itsRadial ->GetType();
    atype=itsAngular->GetType();
    nr=itsNRadial ->get_value_as_int();
    na=itsNAngular->get_value_as_int();
    archive( rtype,atype,nr,na); // serialize things by passing them to the archive
  }
  template<class Archive>  void load(Archive & archive)
  {
    archive( rtype,atype,nr,na); 
    itsRadial  ->set_selected(rtype);
    itsAngular ->set_selected(atype);
    itsNRadial ->set_value(nr);
    itsNAngular->set_value(na);
  }
private:
 
  enumDropDown<qchem::RadialType>* itsRadial;
  enumDropDown<qchem::AngleType>* itsAngular;
  Gtk::SpinButton* itsNRadial;
  Gtk::SpinButton* itsNAngular;

  mutable qchem::RadialType rtype;
  mutable qchem::AngleType  atype;
  mutable guint nr,na;
};

#endif //_MeshFrame_H_