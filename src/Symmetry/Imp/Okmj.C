module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <algorithm>
module qchem.Symmetry.Okmj;

// import qchem.Symmetry.Yl; 
import qchem.Symmetry.AtomEC;
import qchem.Common.Strings; //To get SPDFG string table.
import qchem.stl_io;

using std::cout;
using std::endl;

Omega_k_Sym::Omega_k_Sym(int _kappa) : kappa(_kappa) 
{
    assert(abs(kappa)<10);
};

size_t Omega_k_Sym::SequenceIndex() const //Used for op<
 {
    assert(abs(kappa)<=LMax+1);
    return kappa+LMax;
 }


int Omega_k_Sym::GetDegeneracy() const
{
    return Getj()+0.5; //(2j+1)/2 degeneracy for one spin state.
}


std::ostream& Omega_k_Sym::Write(std::ostream& os) const
{
    int jindex=Getj()-0.5;
    os << SPDFG[GetL()] << j2s[jindex] << " kappa=" << std::setw(2) << kappa << " ";
        
    return os;
}


Omega_kmj_Sym::Omega_kmj_Sym(int _kappa, const std::vector<double>& _mjs) : kappa(_kappa), mjs(_mjs) {};

size_t Omega_kmj_Sym::SequenceIndex() const //Used for op<
 {
    assert(abs(kappa)<=LMax+1);
    double mjmax=*std::max_element(mjs.begin(),mjs.end());
    size_t offset=2*(LMax+1); //End of Omega_k_Sym sequence indexes.
    for (int k1=-(int)LMax-1;k1<kappa;k1++) offset+=2*j(k1)+1; //add up all the degneracies below kappa.
    return mjmax+Getj()+offset;
 }


int Omega_kmj_Sym::GetDegeneracy() const
{
    return mjs.size();
}

std::ostream& Omega_kmj_Sym::Write(std::ostream& os) const
{
    int jindex=Getj()-0.5;
    os << SPDFG[GetL()] << j2s[jindex] << " kappa=" << std::setw(2) << kappa << " mj=" << std::setw(4) << std::setprecision(1) << mjs << " ";
    return os;
}



