// File: HamiltonianFrame.H  GTK frame to show and manage Hamiltonian settings.
#ifndef _HamiltonianFrame_H_
#define _HamiltonianFrame_H_

#include <gtkmm.h>
#include <types.H>

class PlotWindow;

class HamiltonianFrame : public Gtk::Frame
{
  enum htypes {H1E,HF,DFT,D1E,DHF};
public:
  typedef std::shared_ptr<const Cluster> cl_t;
  HamiltonianFrame();
  HamiltonianFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~HamiltonianFrame();

  htypes GetType() const;
  Hamiltonian*  create(const cl_t& cl,const MeshParams* m, const BasisSet* bs) const;
  WaveFunction* create(BasisSet*, ElectronConfiguration*) const; //Decides polarized or un-polarized
  PlotWindow*   create_orbital_pw(BasisSet*,WaveFunction*) const; //Decides polarized or un-polarized
  
  template<class Archive>  void serialize(Archive & archive)
  {
    h_type=GetType();
    is_polarized=itsPolarized->get_active();
    archive( h_type,is_polarized ); // serialize things by passing them to the archive
  }
  void init();

private:
  
  static const std::map<Glib::ustring,htypes> htype_map;
  static std::map<htypes,guint> htype_invmap;
  static htypes find(Glib::ustring);

  const Glib::RefPtr<Gtk::Builder> itsBuilder;
  Gtk::DropDown* itsType;
  Glib::RefPtr<Gtk::StringList> itsTypes; 
  Gtk::CheckButton* itsPolarized;

  mutable htypes h_type;
  mutable bool   is_polarized;
};

#endif //_HamiltonianFrame_H_
