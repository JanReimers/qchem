// File: PolarizedWF.C  Wave function for an unpolarized atom.
module;
#include <cassert>
#include <iostream>
#include <iomanip>
#include <set>
#include "tabulate/table.hpp"

module qchem.WaveFunction.Internal.PolarizedWF;
import qchem.Math;                 // std::round/abs for the fractional-occupation formatting
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
    return new qchem::ChargeDensity::tSpinDensity<T>(GetChargeDensity(Spin::Up),GetChargeDensity(Spin::Down));
}


template <class T> void tPolarizedWF<T>::DisplayEigen() const
{
    Table eigen_table;
    eigen_table.format().multi_byte_characters(true);
    eigen_table.add_row({"Occ/Degen ↑","ϵ↑ (au)","n,Symmetry","Occ/Degen ↓","ϵ↓ (au)","ϵ↑-ϵ↓ (au)"});


    EnergyLevels els_up=this->GetEnergyLevels(Spin::Up), els_dn=this->GetEnergyLevels(Spin::Down);
    // The HIGHEST occupied energy over BOTH channels (see the tUnPolarizedWF twin): a doubly-empty level is
    // dropped only when it sits ABOVE the frontier, where it is one virtual among many.  BELOW the frontier it
    // is a HOLE -- the whole point of looking -- and dropping it made the non-aufbau structure of a MOM-pinned
    // run invisible.  NB the old rule was silently INCONSISTENT between cold and smeared runs: under Fermi
    // smearing no occupation is exactly 0.0, so nothing was ever dropped; at kT=0 the deep empty levels
    // vanished.  Same table, same run, different visibility depending on kT.
    double eHomo=-1e300;
    for (auto elp:this->GetEnergyLevels()) if (elp.second.occ > 0.0) eHomo=std::max(eHomo,elp.second.e);
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
        // ABSENT IS NOT EMPTY (2026-08-10, MnO run 29).  These used to fall back to `el.e` -- the COMBINED
        // level's energy, i.e. the OTHER channel's number -- and to a fabricated occupancy of 0.  The row then
        // read as a level that is occupied in one channel and EMPTY AT THE SAME ENERGY in the other, with
        // ϵ↑−ϵ↓ printing exactly 0.00000000; on MnO that manufactured a spin-up "hole" at −0.0372 below
        // occupied spin-up levels at +0.214, and the campaign's non-aufbau readings came off rows like it.
        // Two different spin Fock matrices do not share an eigenvalue to 8 digits wherever m≠0, so an exact
        // zero in that column was always the tell.  Now an absent level prints "--" and contributes no
        // difference.  (ROOT CAUSE, deeper and still open: `n` indexes a DEGENERATE GROUP, and the grouping
        // differs between channels once ϵ↑≠ϵ↓ -- so pairing rows by (n,sym) is not a sound identity in a
        // polarized run.  It is also why the table skips indices and runs them out of order.  The real fix is
        // to pair by energy/character; this one only stops the display from inventing what it lacks.)
        const int    upDeg=up?up->degen:(dn?dn->degen:1), dnDeg=dn?dn->degen:upDeg;
        if (upOcc==0.0 && dnOcc==0.0 && el.e>eHomo) continue;
        std::ostringstream sym_string,up_occ_string,dn_occ_string;
        sym_string << el.qns.n << *el.qns.sym;
        // Integer occ (gapped insulator) unchanged; FRACTIONAL (Fermi-smeared) shown with decimals, mirroring
        // the unpolarized table.  setprecision(0) alone rounded a smeared 0.996 to "1/1" and 0.004 to "0/1",
        // i.e. it printed a clean integer configuration for a run whose own trace column was flagging partial
        // occupancy every iteration -- the table has to agree with the `m` flag beside it.
        auto occStr=[](std::ostringstream& os, double occ, int deg)
        {
            if (fabs(occ-round(occ)) < 1e-6) os << std::fixed << std::setprecision(0) << occ;
            else                             os << std::fixed << std::setprecision(2) << occ;
            os << "/" << deg;
        };
        if (up) occStr(up_occ_string, upOcc, upDeg); else up_occ_string << "--";
        if (dn) occStr(dn_occ_string, dnOcc, dnDeg); else dn_occ_string << "--";
        size_t l=el.qns.sym->GetPrincipleOffset();
        //! This channel has no such level -> "--", never a number borrowed from the other channel.
        auto eStr=[](const Orbitals::EnergyLevel* lv)
        {
            if (!lv) return std::string("--");
            std::ostringstream os; os << std::fixed << std::setprecision(8) << lv->e; return os.str();
        };
        std::string dEStr="--";
        if (up && dn) { std::ostringstream os; os << std::fixed << std::setprecision(8) << up->e-dn->e; dEStr=os.str(); }

        RowStream rs;
        rs << up_occ_string.str() << eStr(up);
        rs << sym_string.str();
        rs << dn_occ_string.str() << eStr(dn);
        rs << dEStr;
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
template class tPolarizedWF<dcmplx>;

} //namespace
