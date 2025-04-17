// File: SCFFrame.H  GTK frame to show and manage SCF iteration paramater settings.
#ifndef _SCFFrame_H_
#define _SCFFrame_H_

#include <gtkmm.h>
#include <types.H>

class SCFIterationParamsFrame : public Gtk::Frame
{
public:
    SCFIterationParamsFrame() {};
    SCFIterationParamsFrame(BaseObjectType* cobject, const Glib::RefPtr<Gtk::Builder>& refBuilder);
    virtual ~SCFIterationParamsFrame() {};

    SCFIterationParams create() const;
private:
    Gtk::SpinButton*  itsNIter;
    Gtk::Entry*       itsMinDeltaRo;
    Gtk::Entry*       itsRoRelax;
    Gtk::CheckButton* itsVerbose;
};

#endif //_SCFFrame_H_