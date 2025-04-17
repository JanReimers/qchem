// File: LAParamsFrame.cpp  GTK frame to show and manage linear algebra settings.

#include "LAParamsFrame.H"
#include <iostream>


const std::map<Glib::ustring,qchem::Pkg> LAParamsFrame::pkgtype_map=
{
  {"LAPack",qchem::Lapack},
  {"OML"   ,qchem::OML},
};

const std::map<Glib::ustring,qchem::Ortho> LAParamsFrame::orthotype_map=
{
  {"Cholsky",qchem::Cholsky},
  {"SVD"    ,qchem::SVD},
  {"Eigen"  ,qchem::Eigen}, 
};

LAParamsFrame::LAParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder) 
  : Glib::ObjectBase("laparams_frame")
  , Gtk::Frame(cobject)
  , itsPkgType(refBuilder->get_widget<Gtk::DropDown>("la_pkg"))
  , itsOrthoType(refBuilder->get_widget<Gtk::DropDown>("la_ortho"))
  , itsTruncationTolerance(refBuilder->get_widget<Gtk::Entry>("la_trunc"))
  , itsAbsTol(refBuilder->get_widget<Gtk::Entry>("la_abstol"))
{
  std::vector<Glib::ustring> strings;
  for (const auto& [key, _] : pkgtype_map) strings.push_back(key);
  itsPkgTypes = Gtk::StringList::create(strings);
  itsPkgType->set_model(itsPkgTypes);
  itsPkgType->set_selected(0);

  strings.clear();
  for (const auto& [key, _] : orthotype_map) strings.push_back(key);
  itsOrthoTypes = Gtk::StringList::create(strings);
  itsOrthoType->set_model(itsOrthoTypes);
  itsOrthoType->set_selected(0);
}

qchem::Pkg LAParamsFrame::find_pkg(Glib::ustring s)
{
    auto i=pkgtype_map.find(s);
    if (i==pkgtype_map.end())
    {
      std::cerr << "LAParamsFrame::find Unknown package type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

qchem::Ortho LAParamsFrame::find_ortho(Glib::ustring s)
{
    auto i=orthotype_map.find(s);
    if (i==orthotype_map.end())
    {
      std::cerr << "LAParamsFrame::find Unknown ortho type '" << s << "'" << std::endl;
      exit(-1);
    }
    return i->second;
}

LAParams LAParamsFrame::create() const
{
  guint it=itsPkgType->get_selected();
  Glib::ustring pkg_stype=itsPkgTypes->get_string(it);  
  qchem::Pkg  pkg_type=find_pkg(pkg_stype);
  it=itsOrthoType->get_selected();
  Glib::ustring ortho_stype=itsOrthoTypes->get_string(it);  
  qchem::Ortho  ortho_type=find_ortho(ortho_stype);

  double trunctol=Glib::Ascii::strtod(itsTruncationTolerance->get_text());
  double abstol=Glib::Ascii::strtod(itsAbsTol->get_text());
  return {pkg_type,ortho_type,trunctol,abstol};
}

