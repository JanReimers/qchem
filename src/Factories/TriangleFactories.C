// File: TriangleFactories.C

#include "Misc/Triangle/Triangle.H"
#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>


Triangle* Triangle::Factory(std::istream& is)
{
    std::string Name=PeekAtName(is);
    if (Name==typeid(Triangle).name()) return new Triangle;

    std::cout << "Unknown Triangle function type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

