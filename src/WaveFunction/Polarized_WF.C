// File: Polarized_WF.C  Wave function for an unpolarized atom.

#include <cassert>
#include <iostream>
#include <iomanip>
#include "tabulate/table.hpp"

#include "Polarized_WF.H"
#include <ChargeDensity/ChargeDensity.H>
#include <ChargeDensity/Factory.H>
import qchem.Symmetry;

using std::cout;
using std::endl;


Polarized_WF::Polarized_WF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc)
    : Composite_WF(bs,ec,acc) 
{
    MakeIrrep_WFs(Spin::Up);
    MakeIrrep_WFs(Spin::Down);
};

DM_CD* Polarized_WF::GetChargeDensity() const
{
    return PolarizedCD_Factory(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}

WaveFunction::sf_t* Polarized_WF::GetSpinDensity() const
{
    return new SpinDensity(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}

using namespace tabulate;

Color l_colors[]={Color::none,Color::cyan,Color::magenta ,Color::red};
void Polarized_WF::DisplayEigen() const
{
    Table eigen_table;
    eigen_table.format().multi_byte_characters(true);
    eigen_table.add_row({"Occ/Degen ↑","ϵ↑ (au)","n,Symmetry","Occ/Degen ↓","ϵ↓ (au)","ϵ↑-ϵ↓ (au)"});


    EnergyLevels els_up=GetEnergyLevels(Spin::Up), els_dn=GetEnergyLevels(Spin::Down);
    for (auto [eup,up]:els_up)
    {
        if (eup>0.0 || up.occ==0) break;
        Orbital_QNs dnqns(up.qns.n,Spin::Down,up.qns.sym);
        auto dn=els_dn.find(dnqns); 
        std::ostringstream sym_string,up_occ_string,dn_occ_string;
        sym_string << up.qns.n << *up.qns.sym;
        up_occ_string << std::fixed << std::setprecision(0) << up.occ << "/" << up.degen;
        dn_occ_string << std::fixed << std::setprecision(0) << dn.occ << "/" << dn.degen;
        size_t l=up.qns.sym->GetPrincipleOffset();

        RowStream rs;
        rs << up_occ_string.str() << std::fixed << std::setprecision(8) << eup; 
        rs << sym_string.str();
        rs << dn_occ_string.str() << std::fixed << std::setprecision(8) << dn.e;
        rs << std::fixed << up.e-dn.e;
        eigen_table.add_row(rs);
        // Row formating.
        size_t n=eigen_table.size()-1;
        eigen_table[n].format().font_color(l_colors[l]);
        if (dn.occ==0.0)
            for (size_t i:{3,4,5})
            {
                eigen_table[n][i].format().font_style({FontStyle::dark});//.hide_border_top();
                if (n>1) eigen_table[n][i].format().hide_border_top();
            }
        

    }
    // Final table formating.
    size_t N=eigen_table.size();
    for (size_t i=1;i<N-1;i++) eigen_table[i].format().hide_border_bottom();
    for (size_t i=2;i<N;i++) eigen_table[i].format().hide_border_top();
    for (size_t i:{1,4,5}) eigen_table.column(i).format().font_align(FontAlign::right);
    for (size_t i:{0,2,3}) eigen_table.column(i).format().font_align(FontAlign::center);
    cout << eigen_table << endl;
}

