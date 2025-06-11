// File: SCFFrame.cpp  GTK frame to show and manage SCF iteration paramater settings.

#include "SCFFrame.H"
#include <SCFParams.H>

SCFIterationParamsFrame::SCFIterationParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("iteration_frame")
  , Gtk::Frame(cobject)
  , itsNIter(refBuilder->get_widget<Gtk::SpinButton>("SCF_max_iter"))
  , itsMinDeltaRo(refBuilder->get_widget<Gtk::Entry>("SCF_delta_ro"))
  , itsRoRelax(refBuilder->get_widget<Gtk::Entry>("SCF_ro_relax"))
  , itsVerbose(refBuilder->get_widget<Gtk::CheckButton>("SCF_verbose"))
{
  
}

SCFParams SCFIterationParamsFrame::create() const
{
  guint N_iter=itsNIter->get_value_as_int();
  double min_delta_ro=Glib::Ascii::strtod(itsMinDeltaRo->get_text());
  double ro_relax=Glib::Ascii::strtod(itsRoRelax->get_text());
  bool verbose=itsVerbose->get_active();
  return {N_iter,min_delta_ro,ro_relax,0.0,verbose};
}

