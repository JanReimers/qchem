// File: enumDropDown.H  Custom dropdown that works with enums.
#ifndef _enumDropDown_H_
#define _enumDropDown_H_

#include <gtkmm.h>
#include <iostream>
//
//  Isloate all the complexities of dropdowns for enum lists.
//
template <class E> class enumDropDown : public Gtk::DropDown
{
public:
  enumDropDown() {};
  enumDropDown(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder)
  : Glib::ObjectBase(nullptr)
  , Gtk::DropDown(cobject)
  {}

  //  Load the lists of enums and thier corresponding string labels.
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
        std::cerr << "enumDropDown::find Unknown enum type '" << s << "'" << std::endl;
        exit(-1);
    }
    return i->second;
  }
  std::map<Glib::ustring,E> sE_map; //translate string into enum.
  std::map<E,guint> Eu_map; //Translate enum into selector index
  Glib::RefPtr<Gtk::StringList> itsStrings;
};

#endif //_enumDropDown_H_