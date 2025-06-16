// File: json.C  Test out basis set factory using the json header library 

#include "gtest/gtest.h"
#include <BasisSet/BasisSet.H>
#include <BasisSet/Factory.H>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

class jsonTests : public ::testing::Test
{
public:
    jsonTests()
    {
        StreamableObject::SetToPretty();
    }
};

TEST_F(jsonTests,Slater)
{
    json js = {
        {"N", 10},
        {"emin", 0.1},
        {"emax", 5000.0},
    };

    for (size_t Z=1;Z<=92;Z++)
        std::cout << *BasisSetAtom::Factory(BasisSetAtom::Type::Slater,js,Z) << std::endl;

}

TEST_F(jsonTests,Gaussian)
{
    json js = {
        {"N", 9},
        {"emin", 0.1},
        {"emax", 5000.0},
    };

    for (size_t Z=1;Z<=92;Z++)
        std::cout << *BasisSetAtom::Factory(BasisSetAtom::Type::Gaussian,js,Z) << std::endl;
}

TEST_F(jsonTests,BSpline)
{
    json js = {
        {"N", 20},
        {"rmin", 0.1},
        {"rmax", 50.0},
    };

    for (size_t Z=1;Z<=92;Z++)
    {
        std::cout << "Z = " << Z << std::endl;
        std::cout << *BasisSetAtom::Factory(BasisSetAtom::Type::BSpline, js,Z) << std::endl;

    }
}