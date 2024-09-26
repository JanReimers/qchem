// File: BasisSetFactories.C

#include "BasisSetImplementation/PlaneWave/PlaneWaveIE.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  Integral engine factory, reads in name of IntegralEngine derived
//  class and makes a new object using the default constructor.
//

template <class T> IntegralEngine<T>* IntegralEngine<T>::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(        PlaneWaveIE).name()) return new         PlaneWaveIE;

    std::cout << "Unknown integral engine type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

template IntegralEngine<std::complex<double> >* IntegralEngine<std::complex<double> >::Factory(std::istream& is);

