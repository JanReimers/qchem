// File: LAParamsFrame.H  GTK frame to show and manage linear algebra settings.
#ifndef _LAParamsFrame_H_
#define _LAParamsFrame_H_

#include "enumDropDown.H"
#include <gtkmm.h>
import qchem.LAParams;

class LAParamsFrame : public Gtk::Frame
{
public:
  LAParamsFrame() {};
  LAParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~LAParamsFrame() {};

  LAParams create() const;

  template<class Archive>  void save(Archive & archive) const
  {
    ortho=itsOrtho->GetType();
    trunc_tol=Glib::Ascii::strtod(itsTruncationTolerance->get_text());
    archive( ortho, trunc_tol); // serialize things by passing them to the archive
  }
  template<class Archive>  void load(Archive & archive)
  {
    archive( pkg,ortho, trunc_tol, abs_tol); 
    itsPkg->set_selected(pkg);
    itsOrtho->set_selected(ortho);
    itsTruncationTolerance->set_text(Glib::Ascii::dtostr(trunc_tol));
  }

  qchem::Ortho GetOrtho() const;

private:
  enumDropDown<qchem::Ortho>* itsOrtho;
  Gtk::Entry*    itsTruncationTolerance;

  mutable qchem::Ortho ortho;
  mutable double trunc_tol;
};

#endif //_LAParamsFrame_H_