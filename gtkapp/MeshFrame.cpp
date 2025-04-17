// File: MeshFrame.cpp  GTK frame to show and manage DFT integration mesh settings.

#include "MeshFrame.H"
#include "Imp/Mesh/LogRadialMesh.H"
#include "Imp/Mesh/MHLRadialMesh.H"
#include "Imp/Mesh/GaussAngularMesh.H"
#include "Imp/Mesh/GaussLegendreAngularMesh.H"
#include "Imp/Mesh/EulerMaclarenAngularMesh.H"
#include "Imp/Cluster/AtomMesh.H"
#include <iostream>


const std::map<Glib::ustring,qchem::RadialType> MeshFrame::radial_type_map=
{
  {"MHL",qchem::MHL},
  {"Log",qchem::Log},
};
const std::map<Glib::ustring,qchem::AngleType> MeshFrame::angular_type_map=
{
  {"Gauss"         ,qchem::Gauss},
  {"Gauss Legendre",qchem::GaussLegendre},
  {"Euler Mclaren" ,qchem::EulerMclaren},
};

MeshFrame::MeshFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("mesh_frame")
  , Gtk::Frame(cobject)
  , itsRadialType (refBuilder->get_widget<Gtk::DropDown>("mesh_radial"))
  , itsAngularType(refBuilder->get_widget<Gtk::DropDown>("mesh_angular"))
  , itsNRadial    (refBuilder->get_widget<Gtk::SpinButton>("mesh_nr"))
  , itsNAngular   (refBuilder->get_widget<Gtk::SpinButton>("mesh_na"))
{
  std::vector<Glib::ustring> strings;
  for (const auto& [key, _] : radial_type_map) strings.push_back(key);
  itsRadialTypes = Gtk::StringList::create(strings);
  itsRadialType->set_model(itsRadialTypes);
  itsRadialType->set_selected(0);

  strings.clear();
  for (const auto& [key, _] : angular_type_map) strings.push_back(key);
  itsAngularTypes = Gtk::StringList::create(strings);
  itsAngularType->set_model(itsAngularTypes);
  itsAngularType->set_selected(0);
}

qchem::RadialType MeshFrame::find_radial(Glib::ustring s)
{
    auto i=radial_type_map.find(s);
    if (i==radial_type_map.end())
    {
      std::cerr << "MeshFrame::find Unknown radial mesh type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

qchem::AngleType MeshFrame::find_angular(Glib::ustring s)
{
    auto i=angular_type_map.find(s);
    if (i==angular_type_map.end())
    {
      std::cerr << "MeshFrame::find Unknown angular mesh type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

MeshParams MeshFrame::create() const
{
  guint it=itsRadialType->get_selected();
  Glib::ustring r_stype=itsRadialTypes->get_string(it);  
  qchem::RadialType r_type=find_radial(r_stype);
  it=itsAngularType->get_selected();
  Glib::ustring a_stype=itsAngularTypes->get_string(it);  
  qchem::AngleType a_type=find_angular(a_stype);

  guint Nr=itsNRadial->get_value_as_int();
  guint Na=itsNAngular->get_value_as_int();
  int m=2,L=Na;
  double alpha=2.0;

  return {r_type,Nr,m,alpha,a_type,Na,1,1,1};
}

