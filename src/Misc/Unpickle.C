// File: Unpickle.C

#include "BasisSet.H"
#include "Hamiltonian.H"
#include "ChargeDensity.H"
#include "FittedCD.H"
#include "Cluster.H"
#include "Cluster/UnitCell.H"
#include "Mesh/Mesh.H"
#include "WaveFunction.H"
#include <fstream>
#include <iostream>
#include <string>

template <class T> bool UnPickle(T*& pointer, const  char* filep, const char* name)
{
    std::string file(filep);
    bool file_error=true;
    if(file !="")
    {
        std::ifstream in(file.c_str());
        if(!in)
        {
            std::cerr << "Can't open " << name << " file :" << file << std::endl;
            file_error=false;
        }
        else
        {
            pointer = T::Factory(in);
            in >> *pointer;
        }
    }
    return file_error;
}

template bool UnPickle(Hamiltonian        *&,const char*,const char*);
template bool UnPickle(BasisSet           *&,const char*,const char*);
template bool UnPickle(ChargeDensity      *&,const char*,const char*);
template bool UnPickle(FittedCD           *&,const char*,const char*);
template bool UnPickle(Cluster            *&,const char*,const char*);
template bool UnPickle(UnitCell           *&,const char*,const char*);
template bool UnPickle(Mesh               *&,const char*,const char*);
template bool UnPickle(WaveFunction       *&,const char*,const char*);

