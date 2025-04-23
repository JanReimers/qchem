// File: BasisSetFrame.H  GTK frame to show and manage basis set settings.
#ifndef _BasisSetFrame_H_
#define _BasisSetFrame_H_

#include "enumDropDown.H"
#include <gtkmm.h>
#include <types.H>

class BasisSetFrame : public Gtk::Frame
{
public:
  BasisSetFrame();
  BasisSetFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~BasisSetFrame();

  BasisSet* create(size_t LMax) const;

  template<class Archive>  void save(Archive & archive) const
  {
    bstype=itsEnumDD->GetType();
    emin=Glib::Ascii::strtod(itsEmin->get_text());
    emax=Glib::Ascii::strtod(itsEmax->get_text());
    N=itsN->get_value_as_int();
    archive( bstype, emin, emax, N); // serialize things by passing them to the archive
  }
  template<class Archive>  void load(Archive & archive)
  {
    archive(  bstype, emin, emax, N);
    itsEnumDD->set_selected(bstype);
    itsEmin->set_text(Glib::Ascii::dtostr(emin));
    itsEmax->set_text(Glib::Ascii::dtostr(emax));
    itsN->set_value(N);
  }

private:
  enum bstypes {SlaterYl,SlaterYlm,GaussianYl,GaussianYlm};

  enumDropDown<bstypes>* itsEnumDD;
  Gtk::Entry* itsEmin;
  Gtk::Entry* itsEmax;
  Gtk::SpinButton* itsN;

  mutable bstypes bstype;
  mutable double emin,emax;
  mutable guint N;
};


#endif //_BasisSetFrame_H_