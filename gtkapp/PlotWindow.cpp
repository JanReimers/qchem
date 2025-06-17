// File: PlotWindow.cpp Define a window for PLPlot.

#include "PlotWindow.H"
#include <gtkmm-plplot/plotobject2dtext.h>
#include <gtkmm-plplot/plotobject2dline.h>

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
#include <BasisSet/BasisSet.H>
#include <Orbitals/Orbitals.H>
#include <Symmetry/Orbital_QNs.H>
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

void PlotWindow::AddLabels(Glib::ustring x_label, Glib::ustring y_label,Glib::ustring title)
{
  plot.set_axis_title_y(y_label);
  plot.set_axis_title_x(x_label);
  plot.set_plot_title(title);
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

SF_PW::SF_PW(const ScalarFunction<double>* sf,Glib::ustring x_label, Glib::ustring y_label,Glib::ustring title)
{
  std::valarray<double> r(100);
  FillPower(r,0.1,10.0);
  std::valarray<double> y=(*sf)(r);
  y*=r;
  y*=r;
  plot.set_axis_logarithmic_x();
  auto data=Gtk::manage(new Gtk::PLplot::PlotData2D(r,y, Gdk::RGBA("Red"), Gtk::PLplot::LineStyle::CONTINUOUS, 2.0));
  plot.add_data(*data);
  plot.set_axis_title_y(y_label);
  plot.set_axis_title_x(x_label);
  plot.set_plot_title(title);
}

const double dx_sym=1.0;
const double dx_arrow=0.1;
const double dx_spin=dx_arrow+0.001;

#include <WaveFunction/EnergyLevel.H>

//
//  Common construction
//
EnergyLevel_PW::EnergyLevel_PW()
: its_els({
  {Spin::None,EnergyLevels()},
  {Spin::Up  ,EnergyLevels()},
  {Spin::Down,EnergyLevels()}
})
{
  plot.set_axis_title_y("Energy (a.u.)");
  plot.set_axis_title_x("Angular Momentum l");
  plot.set_plot_title("Occupied Orbital Energies");
  plot.hide_legend();
}


//
//  Specific to un-polarized
//
EnergyLevel_PW::EnergyLevel_PW(const WaveFunction* wf)
  : EnergyLevel_PW()
{
  EnergyLevels els=wf->GetEnergyLevels();
  LoadSymmetries(els);
  // std::cout << "Nocc,Nup,Ndn=" << Nocc << " " << Nup << " " << Ndn << std::endl;
  assert(Nup==Ndn);
  DrawArrows(its_els[Spin::None],-0.01,"↑");
  DrawArrows(its_els[Spin::None], 0.01,"↓");
  DrawLevels(its_els[Spin::None], 0.0);
 
}

void EnergyLevel_PW::DrawArrows(const EnergyLevels& els,double _dx_spin,Glib::ustring symbol) 
{
  std::vector<double> x,y;
  auto data=Gtk::manage(new Gtk::PLplot::PlotData2D(x,y, Gdk::RGBA("black"), Gtk::PLplot::LineStyle::NONE));

  for (auto iel:els)
  {
    EnergyLevel el=iel.second;
    size_t occ=el.occ;
    assert(occ!=0);
    double x_left=findx(el)-dx_arrow+_dx_spin;
    double dx=2.0*dx_arrow/(2*occ);
    for (size_t n=0;n<occ;n++)
      data->add_datapoint(x_left+(2*n+1)*dx,el.e);
  }
  data->set_symbol(symbol);
  data->set_symbol_color(Gdk::RGBA("blue"));
  plot.add_data(*data);
}
void EnergyLevel_PW::DrawLevels(const EnergyLevels& els,double _dx_spin)
{
  for (auto iel:els)
  {
    EnergyLevel el=iel.second;
    double x_center=findx(el)+_dx_spin;
    auto line=Gtk::manage(new Gtk::PLplot::PlotObject2DLine(x_center-dx_arrow,el.e,x_center+dx_arrow,el.e));
    plot.add_object(*line);
  }  
}
double EnergyLevel_PW::findx(const EnergyLevel& el) const
{
  Symmetry_Wrap sw(el.qns.sym);
  auto i=its_x_map.find(sw);
  assert(i!=its_x_map.end());
  return i->second;
}
//
//  Assign x values for each symmetry (without spin).
//
void EnergyLevel_PW::LoadSymmetries(const EnergyLevels& els)
{
  double e_scale=fabs(els.begin()->first);
  double x=0.0;
  Nocc=0,Nup=0,Ndn=0;
  StreamableObject::SetToPretty();
  for (auto i:els)
  {
    EnergyLevel el=i.second;
    if (el.occ==0) continue;
    if (el.qns.ms==Spin::None) el.occ/=2; //This akes subsequent arrow drawing code much simpler.
    its_els[el.qns.ms.itsState].insert(el); //Dis-aggregate energy levels by spin state.
    Nocc++;
    if (el.qns.ms==Spin::Up  ) Nup+=el.occ;
    if (el.qns.ms==Spin::Down) Ndn+=el.occ;
    if (el.qns.ms==Spin::None) 
    {
      Nup+=el.occ;
      Ndn+=el.occ;
    }
    
    Symmetry_Wrap sw(el.qns.sym);
    // std::cout << "x,sym = " << x << " " << sw.sym->SequenceIndex() << " " << *sw.sym << std::endl;
    auto is=its_x_map.find(sw);
    if (is==its_x_map.end())
    {
       its_x_map[sw]=x;
       auto text=Gtk::manage(new Gtk::PLplot::PlotObject2DText(el.qns.sym->GetLabel(),x-dx_spin,e_scale*0.1 ));
       plot.add_object(*text);
       x+=dx_sym;
      //  std::cout << "Adding x,sym = " << x << " " << *sw.sym << std::endl;
    }

  }
}

Polarized_EnergyLevel_PW::Polarized_EnergyLevel_PW(const WaveFunction* wf)
: EnergyLevel_PW()
{
  EnergyLevels els=wf->GetEnergyLevels();
  LoadSymmetries(els);
  // std::cout << "Nocc,Nup,Ndn=" << Nocc << " " << Nup << " " << Ndn << std::endl;
  DrawArrows(its_els[Spin::Up  ],-dx_spin,"↑");
  DrawArrows(its_els[Spin::Down], dx_spin,"↓");
  DrawLevels(its_els[Spin::Up  ],-dx_spin);
  DrawLevels(its_els[Spin::Down], dx_spin);


}

#include "oml/vector.h"

template <class T> std::valarray<T> to_valarray(const Vector<T>& v)
{
  std::valarray<T> ret(v.size());
  for (int i:v.arr_indices()) ret[i]=v[i];
  return ret;
}

const std::string ortho_titles[3]={"Cholsky Decomposition Diagonals","Eigen Values","Singlar Values"};
const std::string ortho_ynames[3]={"Cholsky Diagonals","Eigen Values","Singular Values"};


Diagonal_PW::Diagonal_PW(const BasisSet* bs,const WaveFunction* wf, qchem::Ortho ortho, Spin s)
{
  Glib::ustring symbol("•");
  BasisSet::symv_t syms=bs->GetSymmetries();
  int N=syms.size();
  float r=1.0,g=0.0,b=0.0,dr= N==1 ? 1.0 : 1.0/(N-1);
  for (auto sym:syms)
  {
    Irrep_QNs qns(s,sym);
    const TOrbitals<double>* tos=dynamic_cast<const TOrbitals<double>*>(wf->GetOrbitals(qns));
    std::valarray<double> diag=to_valarray(tos->Get_BS_Diagonal());
    auto data=Gtk::manage(new Gtk::PLplot::PlotData2D(diag, Gdk::RGBA("black"), Gtk::PLplot::LineStyle::NONE));
    data->set_symbol(symbol);
    data->set_symbol_color(Gdk::RGBA(r,g,b,1.0));
    data->set_symbol_height_scale_factor(0.8);
    data->set_name(sym->GetLabel()); 
    plot.add_data(*data);
    r-=dr;
    b+=dr;
  }
  plot.set_axis_logarithmic_y();
  plot.set_axis_title_x("Diagonal Index");
  plot.set_axis_title_y(ortho_ynames[ortho]);
  plot.set_plot_title("Basis Set Overlap Matrix "+ortho_titles[ortho]);
}