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
void Orbital_PW::AddLines(const BasisSet* bs, const WaveFunction* wf, Spin s, Glib::ustring symbol)
{
  std::valarray<double> x(100);
  FillPower(x,0.1,10.0);
  BasisSet::symv_t Irreps=bs->GetSymmetries();
  bool use_symbols=symbol!="";
  Glib::ustring spin_symbol="";
  if (s==Spin::Up) spin_symbol="↑";
  if (s==Spin::Down) spin_symbol="↓";

  int line=use_symbols ? Gtk::PLplot::LineStyle::NONE : Gtk::PLplot::LineStyle::CONTINUOUS;
  for (auto sym:bs->GetSymmetries())
  {
    Irrep_QNs qns(s,sym);
    int num_unocc=1; //How many un=occupied orbitals to show?
    const TOrbitals<double>* tos=dynamic_cast<const TOrbitals<double>*>(wf->GetOrbitals(qns));
    int N_to_plot=tos->GetNumOccOrbitals()+num_unocc-1;
    if (N_to_plot==0) N_to_plot=1;
    float r=1.0,g=0.0,b=0.0,dr=1.0/N_to_plot;
    for (auto o=tos->beginT();o!=tos->end();o++)
    {
      double alpha = o->IsOccupied() ? 1.0 : 0.4;
      double thickness = o->IsOccupied() ? 2.0 : 1.0;
      auto data=Gtk::manage(new Gtk::PLplot::PlotData2D(x,(**o)(x), Gdk::RGBA(r,g,b,alpha), (Gtk::PLplot::LineStyle)line, thickness));
      
      if (use_symbols) 
      {
        
        data->set_symbol(symbol);
        data->set_symbol_color(Gdk::RGBA(r,g,b,alpha));
        data->set_symbol_height_scale_factor(0.8);
      }
      data->set_name(o->GetLabel()+spin_symbol); //No to keep out of legend.
      plot.add_data(*data);
      r-=dr;
      b+=dr;
      if (!o->IsOccupied())
      {
        num_unocc--;
        if (num_unocc==0) break;
      }
    }
    if (!use_symbols) line++;
    if (line>8) line=1;
  }
  plot.set_axis_logarithmic_x();
}

UnPolarized_Orbital_PW::UnPolarized_Orbital_PW(const BasisSet* bs, const WaveFunction* wf)
: Orbital_PW()
{
  AddLines(bs,wf,Spin::None);
}


Polarized_Orbital_PW::Polarized_Orbital_PW(const BasisSet* bs, const WaveFunction* wf)
: Orbital_PW()
{
  AddLines(bs,wf,Spin::Up);
  AddLines(bs,wf,Spin::Down,"•");
}

SF_PW::SF_PW(const ScalarFunction<double>* sf)
{
  std::valarray<double> r(100);
  FillPower(r,0.1,10.0);
  std::valarray<double> y=(*sf)(r);
  y*=r;
  y*=r;
  plot.set_axis_logarithmic_x();
  auto data=Gtk::manage(new Gtk::PLplot::PlotData2D(r,y, Gdk::RGBA("Red"), Gtk::PLplot::LineStyle::CONTINUOUS, 2.0));
  plot.add_data(*data);
}
