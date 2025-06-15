// File: MeshFrame.cpp  GTK frame to show and manage DFT integration mesh settings.

#include "MeshFrame.H"
#include "Mesh/LogRadialMesh.H"
#include "Mesh/MHLRadialMesh.H"
#include "Mesh/GaussAngularMesh.H"
#include "Mesh/GaussLegendreAngularMesh.H"
#include "Mesh/EulerMaclarenAngularMesh.H"
#include <iostream>


MeshFrame::MeshFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("mesh_frame")
  , Gtk::Frame(cobject)
  , itsRadial  (Gtk::Builder::get_widget_derived<enumDropDown<qchem::RadialType>>(refBuilder, "mesh_radial"))
  , itsAngular (Gtk::Builder::get_widget_derived<enumDropDown<qchem::AngleType >>(refBuilder, "mesh_angular") )
  , itsNRadial (refBuilder->get_widget<Gtk::SpinButton>("mesh_nr"))
  , itsNAngular(refBuilder->get_widget<Gtk::SpinButton>("mesh_na"))
{
  itsRadial->init({qchem::MHL,qchem::Log},{"MHL","Log"});
  itsAngular->init({qchem::Gauss,qchem::GaussLegendre,qchem::EulerMclaren},{"Gauss","Gauss Legendre","Euler Mclaren"});
}

MeshParams MeshFrame::create() const
{
  rtype=itsRadial ->GetType();
  atype=itsAngular->GetType();
  nr=itsNRadial ->get_value_as_int();
  na=itsNAngular->get_value_as_int();
  size_t m=2;
  double alpha=2.0;

  return {rtype,nr,m,alpha,atype,na,1,1,1};
}

