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
    void AddLabels(Glib::ustring x_label, Glib::ustring y_label,Glib::ustring title);

protected:

    typedef Gtk::PLplot::PlotData2D data_t; 
    Gtk::PLplot::Plot2D plot;
    Gtk::PLplot::Canvas canvas;

    std::vector<std::shared_ptr<data_t>> itsDataSets;
};

#include <types.H>

class Orbital_PW : public PlotWindow
{
public:
    void AddLines(const BasisSet*, const WaveFunction*,Spin s, Glib::ustring symbol="");
};

class UnPolarized_Orbital_PW : public Orbital_PW
{
public:
    UnPolarized_Orbital_PW(const BasisSet*, const WaveFunction*);
    
};

class Polarized_Orbital_PW : public Orbital_PW
{
public:
    Polarized_Orbital_PW(const BasisSet*, const WaveFunction*);

};

// Scalar function plot window
class SF_PW : public PlotWindow 
{
public:
    SF_PW(const ScalarFunction<double>*,Glib::ustring x_label, Glib::ustring y_label,Glib::ustring title);

};

#include <Symmetry/Symmetry.H>
#include <Symmetry/Spin.H>
#include <EnergyLevel.H>
class EnergyLevel_PW : public PlotWindow
{
public:
    EnergyLevel_PW();
    EnergyLevel_PW(const WaveFunction*);

    void LoadSymmetries(const EnergyLevels&);
    void DrawArrows(const EnergyLevels&,double dx_spin,Glib::ustring symbol);
    void DrawLevels(const EnergyLevels&,double dx_spin);
    double findx(const EnergyLevel&) const;


    struct Symmetry_Wrap
    {
        Symmetry_Wrap(const Symmetry* _sym) : sym(_sym) {};
        Symmetry_Wrap(const Symmetry_Wrap& sw) : sym(sw.sym->Clone()) {};
        const Symmetry_Wrap& operator=(const Symmetry_Wrap& sw)
        {
            delete sym;
            sym=sw.sym->Clone();
            return *this;
        }
        friend bool operator<(const Symmetry_Wrap& a, const Symmetry_Wrap& b)
        {
            return a.sym->SequenceIndex() < b.sym->SequenceIndex();
        }
        const Symmetry* sym;
    };
    std::map<Symmetry_Wrap,double> its_x_map;
    std::map<Spin::State,EnergyLevels> its_els;
    size_t Nocc,Nup,Ndn;
};

class Polarized_EnergyLevel_PW : public EnergyLevel_PW
{
public:
    Polarized_EnergyLevel_PW(const WaveFunction*);

};

#include <LASolver/LAParams.H>
class Diagonal_PW : public PlotWindow
{
public:
    Diagonal_PW(const BasisSet*,const WaveFunction*, qchem::Ortho,Spin);

};

#endif //_PlotWindow123_H_