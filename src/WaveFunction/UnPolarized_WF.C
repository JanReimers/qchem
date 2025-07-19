// File: UnPolarized_WF.C  Wave function for an unpolarized atom.

#include <iomanip>
#include <iostream>
#include "tabulate/table.hpp"
#include "UnPolarized_WF.H"
import qchem.Symmetry;
import qchem.Streamable;

UnPolarized_WF::UnPolarized_WF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc)
    : Composite_WF(bs,ec,acc)
{
    MakeIrrep_WFs(Spin::None);
};

using namespace tabulate;

extern Color l_colors[];

void UnPolarized_WF::DisplayEigen() const
{
    Table eigen_table;
    eigen_table.format().multi_byte_characters(true);
    eigen_table.add_row({"Occ/Degen","ϵ (au)","Symmetry"});
       
    for (auto [e,el]:GetEnergyLevels())
    {
       
        if (e>0.0 || el.occ==0) break;
        std::ostringstream sym_string,occ_string;
        sym_string << el.qns.n << *el.qns.sym;
        occ_string << std::fixed << std::setprecision(0) << el.occ << "/" << el.degen;
        size_t l=el.qns.sym->GetPrincipleOffset();

        RowStream rs;
        rs << occ_string.str() << std::fixed << std::setprecision(8) << e; 
        rs << sym_string.str();
        eigen_table.add_row(rs);
        // Row formatting
        size_t n=eigen_table.size()-1;
        eigen_table[n].format().font_color(l_colors[l]);
    }
    // Final table formatting.
    size_t N=eigen_table.size();
    for (size_t i=1;i<N-1;i++) eigen_table[i].format().hide_border_bottom();
    for (size_t i=2;i<N;i++) eigen_table[i].format().hide_border_top();
    for (size_t i:{1}) eigen_table.column(i).format().font_align(FontAlign::right);
    for (size_t i:{0,2}) eigen_table.column(i).format().font_align(FontAlign::center);
    std::cout << eigen_table << std::endl;
}
