// File: PlotWindow.cpp Define a window for PLPlot.

#include "PlotWindow.H"

PlotWindow::PlotWindow()
: plot("X")
, canvas()
{
    // set_default_size(400, 200);
    std::valarray<double> x_va(1000), y_va(1000);
    for (unsigned int i = 0 ; i < 1000 ; i++) {
      x_va[i] = 4*M_PI*i/999;
    }
    y_va = sin(x_va);
    plot.add_data(*Gtk::manage(new Gtk::PLplot::PlotData2D(x_va, y_va, Gdk::RGBA("blue"), Gtk::PLplot::LineStyle::LONG_DASH_LONG_GAP, 5.0)));

    canvas.add_plot(plot);
    const int width = 1024, height = 580;
    canvas.set_vexpand(true);
    canvas.set_hexpand(true);
    canvas.set_focusable(true);
    canvas.set_size_request(-1, height);
    Gtk::AspectFrame geometry(Gtk::Align::CENTER, Gtk::Align::CENTER, float(width)/float(height));
    geometry.set_child(canvas);
    set_child(geometry);
    // set_child(canvas);
}