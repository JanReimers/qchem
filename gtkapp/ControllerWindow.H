// File: ControllerWindow.H  GTK main window for the atom GUI app.
#ifndef _ControllerWindow_H_
#define _ControllerWindow_H_

#include <gtkmm.h>

class AtomFrame;
class BasisSetFrame;
class HamiltonianFrame;
class LAParamsFrame;
class SCFIterationParamsFrame;
class MeshFrame;
class SCFIterator;

class ControllerWindow : public Gtk::Window
{
public:
  ControllerWindow(BaseObjectType* cobject,const Glib::RefPtr<Gtk::Builder>& refBuilder);
  void new_model();

private:
  Gtk::Button* itsStartButton;
  Gtk::Button* itsStepButton;
  Gtk::Button* itsPauseButton;

  AtomFrame*               itsAtom;
  BasisSetFrame*           itsBasisSet;
  HamiltonianFrame*        itsHamiltonian;
  LAParamsFrame*           itsLAParams;
  SCFIterationParamsFrame* itsIterationParams;
  MeshFrame*               itsMeshFrame;
  SCFIterator*             itsSCFIterator;
};

#endif //_ControllerWindow_H_