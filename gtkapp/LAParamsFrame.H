// File: LAParamsFrame.H  GTK frame to show and manage linear algebra settings.
#ifndef _LAParamsFrame_H_
#define _LAParamsFrame_H_

#include "enumDropDown.H"
#include <gtkmm.h>
#include <types.H>
#include <LASolver/LAParams.H>

class LAParamsFrame : public Gtk::Frame
{
public:
  LAParamsFrame() {};
  LAParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~LAParamsFrame() {};

  LAParams create() const;

  template<class Archive>  void save(Archive & archive) const
  {
    pkg=itsPkg->GetType();
    ortho=itsOrtho->GetType();
    trunc_tol=Glib::Ascii::strtod(itsTruncationTolerance->get_text());
    abs_tol=Glib::Ascii::strtod(itsAbsTol->get_text());
    archive( pkg,ortho, trunc_tol, abs_tol); // serialize things by passing them to the archive
  }
  template<class Archive>  void load(Archive & archive)
  {
    archive( pkg,ortho, trunc_tol, abs_tol); 
    itsPkg->set_selected(pkg);
    itsOrtho->set_selected(ortho);
    itsTruncationTolerance->set_text(Glib::Ascii::dtostr(trunc_tol));
    itsAbsTol->set_text(Glib::Ascii::dtostr(abs_tol));
  }

  qchem::Ortho GetOrtho() const;

private:
  enumDropDown<qchem::Pkg  >* itsPkg;
  enumDropDown<qchem::Ortho>* itsOrtho;
  Gtk::Entry*    itsTruncationTolerance;
  Gtk::Entry*    itsAbsTol;

  mutable qchem::Pkg pkg;
  mutable qchem::Ortho ortho;
  mutable double trunc_tol;
  mutable double abs_tol;
};

#endif //_LAParamsFrame_H_