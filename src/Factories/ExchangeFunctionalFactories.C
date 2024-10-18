// File: ExchangeFunctionalFactories.C  Interface for all ExchangeFunctional exchange functions.

#include "Imp/Hamiltonian/SlaterExchange.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

ExchangeFunctional*  ExchangeFunctional::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SlaterExchange).name()) return new  SlaterExchange;

    std::cout << "Unknown  ExchangeFunctional type :" << Name << std::endl;
    std::cout << "Looking for                      :" << typeid(SlaterExchange).name() << std::endl;
    exit(-1);
    return NULL;
}
