// File: HamiltonianFrame.H  GTK frame to show and manage Hamiltonian settings.
#ifndef _HamiltonianFrame_H_
#define _HamiltonianFrame_H_

#include <gtkmm.h>
#include <types.H>

class HamiltonianFrame : public Gtk::Frame
{
public:
  typedef std::shared_ptr<const Cluster> cl_t;
  HamiltonianFrame();
  HamiltonianFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~HamiltonianFrame();

  Hamiltonian* create(const cl_t& cl,const MeshParams* m, const BasisSet* bs) const;

private:
  enum htypes {H1E,HF,DFT,D1E,DHF};
  static const std::map<Glib::ustring,htypes> htype_map;
  static htypes find(Glib::ustring);

  const Glib::RefPtr<Gtk::Builder> itsBuilder;
  Gtk::DropDown* itsType;
  Glib::RefPtr<Gtk::StringList> itsTypes; 
  Gtk::CheckButton* itsPolarized;
};

#endif //_HamiltonianFrame_H_
