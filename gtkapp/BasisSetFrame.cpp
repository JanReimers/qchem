// File: BasisSetFrame.cpp  GTK frame to show and manage basis set settings.

#include "BasisSetFrame.H"
#include "Imp/BasisSet/Atom/l/Slater_BS.H"
#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/ml/Gaussian_BS.H"


BasisSetFrame::BasisSetFrame() : Glib::ObjectBase("basisset_frame") {}

BasisSetFrame::BasisSetFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("basisset_frame")
  , Gtk::Frame(cobject)
  , itsEnumDD(Gtk::Builder::get_widget_derived<enumDropDown<bstypes>>(refBuilder, "basisset_dropdown"))
  , itsEmin(refBuilder->get_widget<Gtk::Entry>("basisset_emin"))
  , itsEmax(refBuilder->get_widget<Gtk::Entry>("basisset_emax"))
  , itsN(refBuilder->get_widget<Gtk::SpinButton>("basisset_N"))
{
  itsEnumDD->init({SlaterYl,SlaterYlm,GaussianYl,GaussianYlm} , {"Slater Yl","Slater Ylm","Gaussian Yl","Gaussian Ylm"} );
}
  
BasisSetFrame::~BasisSetFrame() {};

BasisSet* BasisSetFrame::create(size_t LMax) const
{
  bstype=itsEnumDD->GetType();
  emin=Glib::Ascii::strtod(itsEmin->get_text());
  emax=Glib::Ascii::strtod(itsEmax->get_text());
  N=itsN->get_value_as_int();
  BasisSet* bs=0;
  switch (bstype)
  {
    case SlaterYl : 
      bs=new Atoml::Slater::BasisSet(N,emin,emax,LMax);
      break;
    case SlaterYlm : 
      bs=new Atom_ml::Slater::BasisSet(N,emin,emax,LMax);
      break;
    case GaussianYl : 
      bs= new Atoml::Gaussian::BasisSet(N,emin,emax,LMax);
      break;
    case GaussianYlm : 
      bs= new Atom_ml::Gaussian::BasisSet(N,emin,emax,LMax);
      break;
  } 
  return bs;
}
