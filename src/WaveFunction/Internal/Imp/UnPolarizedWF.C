// File: UnPolarizedWF.C  Wave function for an unpolarized atom.
module;
#include <iomanip>
#include <iostream>
#include "tabulate/table.hpp"

module qchem.WaveFunction.Internal.UnPolarizedWF;
import qchem.SCFAccelerator;
import qchem.Symmetry;
import qchem.Streamable;
import qchem.Strings;

namespace qchem::WaveFunction
{

using namespace tabulate;

template <class T> tUnPolarizedWF<T>::tUnPolarizedWF(const tbs_t<T>* bs,const ElectronConfiguration* ec,tSCFAccelerator<T>* acc)
    : tCompositeWF<T>(bs,ec,acc)
{
    this->MakeIrrepWFs(Spin::None);
};



template <class T> void tUnPolarizedWF<T>::DisplayEigen() const
{
    Table eigen_table;
    eigen_table.format().multi_byte_characters(true);
    eigen_table.add_row({"Occ/Degen","ϵ (au)","Symmetry"});
       
    for (auto [e,el]:this->GetEnergyLevels())
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

template class tUnPolarizedWF<double>;
template class tUnPolarizedWF<dcmplx>;

} //namespace
