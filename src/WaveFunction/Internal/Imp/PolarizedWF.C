// File: PolarizedWF.C  Wave function for an unpolarized atom.
module;
#include <cassert>
#include <iostream>
#include <iomanip>
#include <set>
#include "tabulate/table.hpp"

module qchem.WaveFunction.Internal.PolarizedWF;
import qchem.SCFAccelerator;
import qchem.ChargeDensity.Factory;
import qchem.Streamable;
import qchem.Common.Strings;

namespace qchem::WaveFunction
{

using namespace tabulate;

using std::cout;
using std::endl;


PolarizedWF::PolarizedWF(const bs_t* bs,const ElectronConfiguration* ec,SCFAccelerator* acc)
    : CompositeWF(bs,ec,acc) 
{
    MakeIrrepWFs(Spin::Up);
    MakeIrrepWFs(Spin::Down);
};

DM_CD* PolarizedWF::GetChargeDensity() const
{
    return PolarizedCD_Factory(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}

WaveFunction::sf_t* PolarizedWF::GetSpinDensity() const
{
    return new qchem::ChargeDensity::SpinDensity(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}


void PolarizedWF::DisplayEigen() const
{
    Table eigen_table;
    eigen_table.format().multi_byte_characters(true);
    eigen_table.add_row({"Occ/Degen ↑","ϵ↑ (au)","n,Symmetry","Occ/Degen ↓","ϵ↓ (au)","ϵ↑-ϵ↓ (au)"});


    EnergyLevels els_up=GetEnergyLevels(Spin::Up), els_dn=GetEnergyLevels(Spin::Down);
    std::set<Orbital_QNs> alreadyGotIt;
    for (auto elp:GetEnergyLevels())
    {
        const Orbitals::EnergyLevel& el=elp.second;
        Orbital_QNs upqns(el.qns.n,Spin::Up  ,el.qns.sym);
        Orbital_QNs dnqns(el.qns.n,Spin::Down,el.qns.sym);
        auto up=els_up.find(upqns); 
        auto dn=els_dn.find(dnqns); 
        if (alreadyGotIt.find(upqns)!=alreadyGotIt.end())
        {
            assert(alreadyGotIt.find(dnqns)!=alreadyGotIt.end());
            continue;
        }
        alreadyGotIt.insert(upqns);
        alreadyGotIt.insert(dnqns);
        if (up.occ==0 && dn.occ==0) continue;
        std::ostringstream sym_string,up_occ_string,dn_occ_string;
        sym_string << up.qns.n << *up.qns.sym;
        up_occ_string << std::fixed << std::setprecision(0) << up.occ << "/" << up.degen;
        dn_occ_string << std::fixed << std::setprecision(0) << dn.occ << "/" << dn.degen;
        size_t l=up.qns.sym->GetPrincipleOffset();

        RowStream rs;
        rs << up_occ_string.str() << std::fixed << std::setprecision(8) << up.e; 
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

} //namespace
