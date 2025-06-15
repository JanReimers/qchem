// File: json.C  Test out basis set factory using the json header library 

#include "gtest/gtest.h"
#include "Symmetry/Atom_EC.H"
#include <BasisSet.H>
#include "Imp/BasisSet/Atom/l/Slater_BS.H"
#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/ml/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/l/BSpline_BS.H"
#include "Imp/BasisSet/Atom/ml/BSpline_BS.H"
#include <cassert>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

enum class AtomBasisType {Slater,Gaussian,BSpline};
enum class AtomAngularBasisType {Yl,Ylm};

namespace BasisSetAtom
{
    BasisSet* Factory(const nlohmann::json& js,const ElectronConfiguration& ec)
    {
        AtomBasisType type=js["type"].template get<AtomBasisType>();
        size_t N=js["N"];
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax(),g=2*LMax+1;
        ml_Breakdown bd=aec.GetBreadown(LMax);
        AtomAngularBasisType atype = (bd.ml_paired.size()==g || bd.ml_unpaired.size()==g) ? AtomAngularBasisType::Yl : AtomAngularBasisType::Ylm;
        BasisSet* bs=0;
        switch (type)
        {
        case AtomBasisType::Slater:
        {
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            switch (atype)
            {
            case AtomAngularBasisType::Yl:
                bs=new Atoml::Slater::BasisSet(N,emin,emax,LMax);
                break;
            
            case AtomAngularBasisType::Ylm:
                bs=new Atom_ml::Slater::BasisSet(N,emin,emax,ec);
            }
            break;
        }
        case AtomBasisType::Gaussian:
        {
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            switch (atype)
            {
            case AtomAngularBasisType::Yl:
                bs=new Atoml::Gaussian::BasisSet(N,emin,emax,LMax);
                break;
            
            case AtomAngularBasisType::Ylm:
                bs=new Atom_ml::Gaussian::BasisSet(N,emin,emax,ec);
            }
            break;
        }
        case AtomBasisType::BSpline:
        {
            double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
            switch (atype)
            {
            case AtomAngularBasisType::Yl:
                bs=new Atoml::BSpline::BasisSet<6>(N,rmin,rmax,LMax);
                break;
            
            case AtomAngularBasisType::Ylm:
                bs=new Atom_ml::BSpline::BasisSet<6>(N,rmin,rmax,ec);
            }
            break;   
        } 
        }
    
    assert(bs);
    return bs;
    }
}

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
        {"type", AtomBasisType::Slater},
        {"N", 10},
        {"emin", 0.1},
        {"emax", 5000.0},
    };

    for (size_t Z=1;Z<=92;Z++)
    {
        Atom_EC ec(Z);
        std::cout << *BasisSetAtom::Factory(js,ec) << std::endl;

    }
}

TEST_F(jsonTests,Gaussian)
{
    json js = {
        {"type", AtomBasisType::Gaussian},
        {"N", 9},
        {"emin", 0.1},
        {"emax", 5000.0},
    };

    for (size_t Z=1;Z<=92;Z++)
    {
        Atom_EC ec(Z);
        std::cout << *BasisSetAtom::Factory(js,ec) << std::endl;

    }
}

TEST_F(jsonTests,BSpline)
{
    json js = {
        {"type", AtomBasisType::BSpline},
        {"N", 20},
        {"rmin", 0.1},
        {"rmax", 50.0},
    };

    for (size_t Z=1;Z<=92;Z++)
    {
        Atom_EC ec(Z);
        std::cout << *BasisSetAtom::Factory(js,ec) << std::endl;

    }
}