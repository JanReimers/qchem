// File Slater_m/BasisSet.H

#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
#include "Imp/BasisSet/Atom/ml/Slater_IBS.H"
#include "Imp/BasisSet/Atom/radial/Slater/ExponentScaler.H"
#include "Imp/BasisSet/Atom/EC.H"
#include "Imp/Containers/stl_io.h"
#include <iostream>

using std::cout;
using std::endl;

namespace Atom_ml
{
namespace Slater
{


BasisSet::BasisSet(size_t N, double emin, double emax, size_t LMax)
{
    ::Slater::ExponentScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        for (int m=-L;m<=(int)L;m++)
            Insert(new Orbital_IBS(this,ss.Get_es(L),L,m));            

}

BasisSet::BasisSet(size_t N, double emin, double emax, const ElectronConfiguration& ec)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    int LMax=aec.GetLMax();
    ::Slater::ExponentScaler ss(N,emin,emax,LMax);
    for (int L=0;L<=LMax;L++)
    {
        if (aec.GetNval(L)==0)
        {
            std::vector<int> mls;
            for (int ml=-L;ml<=L;ml++) mls.push_back(ml);
            Insert(new Orbital_IBS(this,ss.Get_es(L),L,mls));   
        }
        else
        {
            size_t g=(2*L+1); //degenracy
            size_t Nunp=aec.GetNUnapired();
            size_t Nl=aec.GetNval(L);
            assert((Nl-Nunp)%2==0);  //This fails Cr Z=24 which one 3s and 5 3d unpaired electrons.
            assert(Nl>=Nunp);
            size_t Npaired= (Nl-Nunp)/2;
            assert(g>=Npaired+Nunp);
            size_t Nempty=g-Npaired-Nunp;
            std::vector<int> ml_paired,ml_unpaired,ml_unoccupied;
            int ml=-(int)L;
            for (size_t i=0;i<Npaired;i++) ml_paired    .push_back(ml++);
            for (size_t i=0;i<Nunp   ;i++) ml_unpaired  .push_back(ml++);
            for (size_t i=0;i<Nempty ;i++) ml_unoccupied.push_back(ml++);
            cout << "ml_paired    =" << ml_paired << endl;
            cout << "ml_unpaired  =" << ml_unpaired << endl;
            cout << "ml_unoccupied=" << ml_unoccupied << endl;
            assert(ml==(int)(L+1));
            if (ml_paired.size()>0)   
                Insert(new Orbital_IBS(this,ss.Get_es(L),L,ml_paired));            
            if (ml_unpaired.size()>0)   
                Insert(new Orbital_IBS(this,ss.Get_es(L),L,ml_unpaired));            
            if (ml_unoccupied.size()>0)   
                Insert(new Orbital_IBS(this,ss.Get_es(L),L,ml_unoccupied));            

        }
    }
    //
    //  TODO move most of this logic into Atom_EC
    //
    //The valvance electrons are not nessexcarily at LMax,  i.e. sodium, valance is 2s.
    
}

}} //namespace
