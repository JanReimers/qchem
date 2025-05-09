// File: LAParamsFrame.cpp  GTK frame to show and manage linear algebra settings.

#include "LAParamsFrame.H"
#include <iostream>

LAParamsFrame::LAParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("laparams_frame")
  , Gtk::Frame(cobject)
  , itsPkg(Gtk::Builder::get_widget_derived<enumDropDown<qchem::Pkg>>(refBuilder, "la_pkg"))
  , itsOrtho(Gtk::Builder::get_widget_derived<enumDropDown<qchem::Ortho>>(refBuilder, "la_ortho"))
  , itsTruncationTolerance(refBuilder->get_widget<Gtk::Entry>("la_trunc"))
  , itsAbsTol(refBuilder->get_widget<Gtk::Entry>("la_abstol"))
{
  itsPkg->init({qchem::OML,qchem::Lapack},{"OML","Lapack"});
  itsOrtho->init({qchem::Cholsky,qchem::Eigen,qchem::SVD},{"Cholsky","Eigen","SVD"});
}


LAParams LAParamsFrame::create() const
{
  pkg=itsPkg->GetType();
  ortho=itsOrtho->GetType();
  trunc_tol=Glib::Ascii::strtod(itsTruncationTolerance->get_text());
  abs_tol=Glib::Ascii::strtod(itsAbsTol->get_text());
  return {pkg,ortho,trunc_tol,abs_tol};
}

qchem::Ortho LAParamsFrame::GetOrtho() const
{
  ortho=itsOrtho->GetType();
  return ortho;
}
