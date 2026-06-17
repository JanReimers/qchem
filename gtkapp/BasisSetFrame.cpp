// File: BasisSetFrame.cpp  GTK frame to show and manage basis set settings.

#include "BasisSetFrame.H"

import qchem.BasisSet.Atom.Factory;

using namespace BasisSet::Atom;

BasisSetFrame::BasisSetFrame() : Glib::ObjectBase("basisset_frame") {}

BasisSetFrame::BasisSetFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("basisset_frame")
  , Gtk::Frame(cobject)
  , itsEnumDD(Gtk::Builder::get_widget_derived<enumDropDown<bstypes>>(refBuilder, "basisset_dropdown"))
  , itsEmin(refBuilder->get_widget<Gtk::Entry>("basisset_emin"))
  , itsEmax(refBuilder->get_widget<Gtk::Entry>("basisset_emax"))
  , itsN(refBuilder->get_widget<Gtk::SpinButton>("basisset_N"))
{
  itsEnumDD->init({Slater,Gaussian,BSpline6,BSpliner6,Gaussian_RKB,Slater_RKB} , {"Slater","Gaussian","BSpline6","BSpliner6","Gaussian_RKB","Slater_RKB"} );
}
  
BasisSetFrame::~BasisSetFrame() {};

BasisSet* BasisSetFrame::create(size_t LMax) const
{
  bstype=itsEnumDD->GetType();
  emin=Glib::Ascii::strtod(itsEmin->get_text());
  emax=Glib::Ascii::strtod(itsEmax->get_text());
  N=itsN->get_value_as_int();
  size_t Z=2;
  nlohmann::json js={{"type",type},{"N",N},{"emin",emin},{"emax",emax}}
  return Factory(js,Z);
}
