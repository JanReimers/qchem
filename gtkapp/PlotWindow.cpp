// File: PlotWindow.cpp Define a window for PLPlot.

#include "PlotWindow.H"

PlotWindow::PlotWindow()
: plot("X")
, canvas()
{
    // set_default_size(400, 200);
    // std::valarray<double> x_va(1000), y_va(1000);
    // for (unsigned int i = 0 ; i < 1000 ; i++) {
    //   x_va[i] = 4*M_PI*i/999;
    // }
    // y_va = sin(x_va);
    // plot.add_data(*Gtk::manage(new Gtk::PLplot::PlotData2D(x_va, y_va, Gdk::RGBA("blue"), Gtk::PLplot::LineStyle::LONG_DASH_LONG_GAP, 5.0)));

    canvas.add_plot(plot);
    const int width = 1024, height = 580;
    canvas.set_vexpand(true);
    canvas.set_hexpand(true);
    canvas.set_focusable(true);
    canvas.set_size_request(width, -1);
    Gtk::AspectFrame geometry(Gtk::Align::CENTER, Gtk::Align::CENTER, float(width)/float(height));
    geometry.set_child(canvas);
    append(geometry);
    // set_child(canvas);
}

#include <valarray>
template <class T> void FillPower(std::valarray<T>& arr,T start, T stop)
{
  double del=0.5*(start+stop); //n=1 case
  if (arr.size()>1)
    del=(std::log(stop/start))/(double)(arr.size()-1);
  auto i=std::begin(arr);
  for (int n=0;i!=std::end(arr);i++,n++) *i=T(start*std::exp(n*del));
}

#include <WaveFunction.H>
#include <BasisSet.H>
#include <Orbital.H>
#include <Orbital_QNs.H>

UnPolarized_Orbital_PW::UnPolarized_Orbital_PW(const BasisSet* bs, const WaveFunction* wf)
: PlotWindow()
{
  std::valarray<double> x(100);
  FillPower(x,0.1,10.0);
  for (auto sym:bs->GetSymmetries())
  {
    Irrep_QNs qns(Spin::None,sym);
    const TOrbitals<double>* tos=dynamic_cast<const TOrbitals<double>*>(wf->GetOrbitals(qns));
    for (auto o=tos->beginT();o!=tos->end();o++)
    {
      plot.add_data(*Gtk::manage(new Gtk::PLplot::PlotData2D(x,(**o)(x), Gdk::RGBA("blue"), Gtk::PLplot::LineStyle::LONG_DASH_LONG_GAP, 1.0)));
    }
  }
  plot.set_axis_logarithmic_x();
}

Polarized_Orbital_PW::Polarized_Orbital_PW(const BasisSet* bs, const WaveFunction* wf)
: PlotWindow()
{
  std::valarray<double> x(100);
  FillPower(x,0.01,40.0);
  for (auto sym:bs->GetSymmetries())
  {
    Irrep_QNs qns(Spin::Up,sym);
    const TOrbitals<double>* tos=dynamic_cast<const TOrbitals<double>*>(wf->GetOrbitals(qns));
    for (auto o=tos->beginT();o!=tos->end();o++)
    {
      plot.add_data(*Gtk::manage(new Gtk::PLplot::PlotData2D(x,(**o)(x), Gdk::RGBA("blue"), Gtk::PLplot::LineStyle::LONG_DASH_LONG_GAP, 1.0)));
    }
  }
}