// File: HamiltonianFrame.H  GTK frame to show and manage Hamiltonian settings.
#ifndef _HamiltonianFrame_H_
#define _HamiltonianFrame_H_

#include <gtkmm.h>
#include <types.H>


template <class E> class enumDropDown : public Gtk::DropDown
{
public:
  enumDropDown() {};
  enumDropDown(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder)
  : Glib::ObjectBase(nullptr)
  , Gtk::DropDown(cobject)
  {}

  void init(std::initializer_list<E> enums,std::vector<Glib::ustring> labels)
  {
    auto s=labels.begin();
    for (auto e:enums) sE_map.insert({*s++,e});
    itsStrings=Gtk::StringList::create(labels);
    for (guint i=0;i<itsStrings->get_n_items();i++)
    {
      E e=find(itsStrings->get_string(i));
      Eu_map[e]=i;
    }
    set_model(itsStrings);
    set_selected(*enums.begin());
  }

  void set_selected(const E& e)
  {
    guint i=Eu_map[e];
    Gtk::DropDown::set_selected(i);
  }

  E GetType() const
  {
    guint it=get_selected();
    Glib::ustring s=itsStrings->get_string(it);  
    return find(s);
  }

private:
  E find(Glib::ustring s) const
  {
    auto i=sE_map.find(s);
    if (i==sE_map.end())
    {
        // std::cerr << "HamiltonianFrame::find Unknown Hamiltonian type '" << s << "'" << std::endl;
        exit(-1);
    }
    return i->second;
  }
  std::map<Glib::ustring,E> sE_map; //translate string into enum.
  std::map<E,guint> Eu_map; //Translate enum into selector index
  Glib::RefPtr<Gtk::StringList> itsStrings;
};

class PlotWindow;

class HamiltonianFrame : public Gtk::Frame
{
  enum htypes {H1E,HF,DFT,D1E,DHF};
public:
  typedef std::shared_ptr<const Cluster> cl_t;
  HamiltonianFrame();
  HamiltonianFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
  virtual ~HamiltonianFrame();

  // htypes GetType() const;
  Hamiltonian*  create(const cl_t& cl,const MeshParams* m, const BasisSet* bs) const;
  WaveFunction* create(BasisSet*, ElectronConfiguration*) const; //Decides polarized or un-polarized
  PlotWindow*   create_orbital_pw(BasisSet*,WaveFunction*) const; //Decides polarized or un-polarized
  
  template<class Archive>  void serialize(Archive & archive)
  {
    h_type=itsEnumDD->GetType();
    is_polarized=itsPolarized->get_active();
    archive( h_type,is_polarized ); // serialize things by passing them to the archive
  }
  void init();

private:
  
  static const std::map<Glib::ustring,htypes> htype_map;
  static std::map<htypes,guint> htype_invmap;
  static htypes find(Glib::ustring);

  const Glib::RefPtr<Gtk::Builder> itsBuilder;
  // Gtk::DropDown* itsType;
  // Glib::RefPtr<Gtk::StringList> itsTypes; 
  Gtk::CheckButton* itsPolarized;
  enumDropDown<htypes>* itsEnumDD;

  mutable htypes h_type;
  mutable bool   is_polarized;
};

#endif //_HamiltonianFrame_H_
