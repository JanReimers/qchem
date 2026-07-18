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
import qchem.Strings;

namespace qchem::WaveFunction
{

using namespace tabulate;

using std::cout;
using std::endl;


template <class T> tPolarizedWF<T>::tPolarizedWF(const tbs_t<T>* bs,const ElectronConfiguration* ec,tSCFAccelerator<T>* acc,
                                                 qchem::Ortho basisOrtho, double basisOrthoTol)
    : tCompositeWF<T>(bs,ec,acc,basisOrtho,basisOrthoTol)
{
    this->MakeIrrepWFs(Spin::Up);
    this->MakeIrrepWFs(Spin::Down);
};

template <class T> tDM_CD<T>* tPolarizedWF<T>::GetChargeDensity() const
{
    return PolarizedCD_Factory(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}

template <class T> typename tPolarizedWF<T>::sf_t* tPolarizedWF<T>::GetSpinDensity() const
{
    return new qchem::ChargeDensity::SpinDensity(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}


template <class T> void tPolarizedWF<T>::DisplayEigen() const
{
    Table eigen_table;
    eigen_table.format().multi_byte_characters(true);
    eigen_table.add_row({"Occ/Degen ↑","ϵ↑ (au)","n,Symmetry","Occ/Degen ↓","ϵ↓ (au)","ϵ↑-ϵ↓ (au)"});


    EnergyLevels els_up=this->GetEnergyLevels(Spin::Up), els_dn=this->GetEnergyLevels(Spin::Down);
    std::set<Orbital_QNs> alreadyGotIt;
    for (auto elp:this->GetEnergyLevels())
    {
        const Orbitals::EnergyLevel& el=elp.second;
        Orbital_QNs upqns(el.qns.n,Spin::Up  ,el.qns.sym);
        Orbital_QNs dnqns(el.qns.n,Spin::Down,el.qns.sym);
        if (alreadyGotIt.find(upqns)!=alreadyGotIt.end()) continue;
        alreadyGotIt.insert(upqns);
        alreadyGotIt.insert(dnqns);
        // A combined level need NOT exist in both spin channels (open shell -- e.g. an O2 triplet level that
        // is occupied/present in ↑ but not ↓).  Guard both (find()==UB on a miss in Release) and take the
        // label + l from el.qns, whose sym is always valid (this fixes the O2-HF-triplet SEGV).
        const Orbitals::EnergyLevel* up=els_up.FindOrNull(upqns);
        const Orbitals::EnergyLevel* dn=els_dn.FindOrNull(dnqns);
        const double upOcc=up?up->occ:0.0, dnOcc=dn?dn->occ:0.0;
        const double upE  =up?up->e:el.e , dnE  =dn?dn->e:el.e;
        const int    upDeg=up?up->degen:(dn?dn->degen:1), dnDeg=dn?dn->degen:upDeg;
        if (upOcc==0.0 && dnOcc==0.0) continue;
        std::ostringstream sym_string,up_occ_string,dn_occ_string;
        sym_string << el.qns.n << *el.qns.sym;
        up_occ_string << std::fixed << std::setprecision(0) << upOcc << "/" << upDeg;
        dn_occ_string << std::fixed << std::setprecision(0) << dnOcc << "/" << dnDeg;
        size_t l=el.qns.sym->GetPrincipleOffset();

        RowStream rs;
        rs << up_occ_string.str() << std::fixed << std::setprecision(8) << upE;
        rs << sym_string.str();
        rs << dn_occ_string.str() << std::fixed << std::setprecision(8) << dnE;
        rs << std::fixed << upE-dnE;
        eigen_table.add_row(rs);
        // Row formating.
        size_t n=eigen_table.size()-1;
        eigen_table[n].format().font_color(l_colors[l]);
        if (dnOcc==0.0)
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

template class tPolarizedWF<double>;

} //namespace
