// File: json.C  Test out basis set factory using the json header library 

#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
import qchem.Factory;
import qchem.BasisSet;
import qchem.Streamable;
using json = nlohmann::json;

class jsonTests : public ::testing::Test
{
public:
    jsonTests()
    {
    }
};

TEST_F(jsonTests,Slater)
{
    json js = {
        {"N", 10},
        {"emin", 0.1},
        {"emax", 5000.0},
        {"type", BasisSetAtomFactory::Type::Slater},
    };

    for (size_t Z=1;Z<=92;Z++)
        BasisSetAtomFactory::Factory(js,Z);

}

TEST_F(jsonTests,Gaussian)
{
    json js = {
        {"N", 9},
        {"emin", 0.1},
        {"emax", 5000.0},
        {"type", BasisSetAtomFactory::Type::Gaussian},
    };

    for (size_t Z=1;Z<=92;Z++)
        BasisSetAtomFactory::Factory(js,Z);

}

TEST_F(jsonTests,BSpline)
{
    json js = {
        {"N", 20},
        {"rmin", 0.1},
        {"rmax", 50.0},
        {"type", BasisSetAtomFactory::Type::BSpline6},
    };

    for (size_t Z=1;Z<=92;Z++)
        BasisSetAtomFactory::Factory(js,Z);
}