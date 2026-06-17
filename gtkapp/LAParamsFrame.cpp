// File: LAParamsFrame.cpp  GTK frame to show and manage linear algebra settings.

#include "LAParamsFrame.H"
#include <iostream>

LAParamsFrame::LAParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("laparams_frame")
  , Gtk::Frame(cobject)
  , itsOrtho(Gtk::Builder::get_widget_derived<enumDropDown<qchem::Ortho>>(refBuilder, "la_ortho"))
  , itsTruncationTolerance(refBuilder->get_widget<Gtk::Entry>("la_trunc"))
{
  itsOrtho->init({qchem::Cholsky,qchem::Eigen,qchem::SVD},{"Cholsky","Eigen","SVD"});
}


LAParams LAParamsFrame::create() const
{
  ortho=itsOrtho->GetType();
  trunc_tol=Glib::Ascii::strtod(itsTruncationTolerance->get_text());
  return {ortho,trunc_tol};
}

qchem::Ortho LAParamsFrame::GetOrtho() const
{
  ortho=itsOrtho->GetType();
  return ortho;
}
