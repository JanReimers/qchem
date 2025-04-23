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

    template<class Archive>  void save(Archive & archive) const
    {
        niter=itsNIter->get_value_as_int();
        min_dro=Glib::Ascii::strtod(itsMinDeltaRo->get_text());
        ro_relax=Glib::Ascii::strtod(itsRoRelax->get_text());
        verbose=itsVerbose->get_active();
        archive( niter, min_dro, ro_relax,verbose); // serialize things by passing them to the archive
    }
    template<class Archive>  void load(Archive & archive)
    {
        archive( niter, min_dro, ro_relax,verbose);
        itsNIter->set_value(niter);
        itsMinDeltaRo->set_text(Glib::Ascii::dtostr(min_dro));
        itsRoRelax->set_text(Glib::Ascii::dtostr(ro_relax));
        itsVerbose->set_active(verbose);
    }

private:
    Gtk::SpinButton*  itsNIter;
    Gtk::Entry*       itsMinDeltaRo;
    Gtk::Entry*       itsRoRelax;
    Gtk::CheckButton* itsVerbose;

    mutable guint niter;
    mutable double min_dro,ro_relax;
    mutable bool verbose;
};

#endif //_SCFFrame_H_