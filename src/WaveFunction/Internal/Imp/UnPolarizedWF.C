// File: UnPolarizedWF.C  Wave function for an unpolarized atom.
module;
#include <iomanip>
#include <iostream>
#include <cmath>
#include "tabulate/table.hpp"

module qchem.WaveFunction.Internal.UnPolarizedWF;
import qchem.SCFAccelerator;
import qchem.Symmetry;
import qchem.Streamable;
import qchem.Strings;

namespace qchem::WaveFunction
{

using namespace tabulate;

template <class T> tUnPolarizedWF<T>::tUnPolarizedWF(const tbs_t<T>* bs,const ElectronConfiguration* ec,tSCFAccelerator<T>* acc,
                                                     qchem::Ortho basisOrtho, double basisOrthoTol)
    : tCompositeWF<T>(bs,ec,acc,basisOrtho,basisOrthoTol)
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
        // Stop past the frontier by OCCUPATION, not energy sign.  The old `e>0.0` cutoff is a MOLECULAR idiom
        // (bound states sit below the vacuum level at 0); in a SOLID the energy zero is arbitrary (the PP
        // G=0/alignment convention), so the Fermi level -- and every occupied level -- can be POSITIVE (a
        // metal: this hid all but the one negative-energy Γ level).  Occupations are monotonic in energy under
        // one μ, so the first negligibly-occupied level ends the occupied set for atoms and metals alike.
        if (el.occ < 1e-6) break;
        std::ostringstream sym_string,occ_string;
        sym_string << el.qns.n << *el.qns.sym;
        // Integer occ (atoms / gapped insulators) unchanged; fractional (Fermi-smeared metal) shown with
        // decimals so a partially-filled band is honest instead of rounding to an integer.
        if (std::abs(el.occ - std::round(el.occ)) < 1e-6)
            occ_string << std::fixed << std::setprecision(0) << el.occ << "/" << el.degen;
        else
            occ_string << std::fixed << std::setprecision(2) << el.occ << "/" << el.degen;
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
