// File: PlotWindow.H Define a window for PLPlot.
#ifndef _PlotWindow123_H_
#define _PlotWindow123_H_

#include <gtkmm.h>
#include "gtkmm-plplot.h"
#include <vector>
#include <memory>

class PlotWindow : public Gtk::Box
{
public:
    PlotWindow();
    virtual ~PlotWindow() {};

protected:

    typedef Gtk::PLplot::PlotData2D data_t; 
    Gtk::PLplot::Plot2D plot;
    Gtk::PLplot::Canvas canvas;

    std::vector<std::shared_ptr<data_t>> itsDataSets;
};

#include <types.H>

class UnPolarized_Orbital_PW : public PlotWindow
{
public:
    UnPolarized_Orbital_PW(const BasisSet*, const WaveFunction*);
};

class Polarized_Orbital_PW : public PlotWindow
{
public:
    Polarized_Orbital_PW(const BasisSet*, const WaveFunction*);

};


#endif //_PlotWindow123_H_