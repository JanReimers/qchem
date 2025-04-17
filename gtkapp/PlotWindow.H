// File: PlotWindow.H Define a window for PLPlot.
#ifndef _PlotWindow123_H_
#define _PlotWindow123_H_

#include <gtkmm.h>
#include "gtkmm-plplot.h"
#include <vector>
#include <memory>

class PlotWindow : public Gtk::Frame
{
public:
    PlotWindow();
    virtual ~PlotWindow() {};

private:

    typedef Gtk::PLplot::PlotData2D data_t; 
    Gtk::PLplot::Plot2D plot;
    Gtk::PLplot::Canvas canvas;

    std::vector<std::shared_ptr<data_t>> itsDataSets;
};


#endif //_PlotWindow123_H_