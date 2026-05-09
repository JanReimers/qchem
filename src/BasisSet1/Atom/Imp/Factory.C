// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.
module; 
#include <cassert>
#include <nlohmann/json.hpp>

module qchem.BasisSet1.Atom.Factory;
import qchem.Symmetry.AtomEC;
// import qchem.BasisSet.Atom.Slater.RKB.BS;
// import qchem.BasisSet.Atom.Gaussian.NR.BS;
// import qchem.BasisSet.Atom.Slater.NR.BS;
// import qchem.BasisSet.Atom.BSpline.NR.BS;

import qchem.BasisSet1.Atom.BasisSet;
import BasisSet.Atom.BSpline.NR.BS_Evaluator;
import BasisSet.Atom.Gaussian_BS;
import BasisSet.Atom.Slater.NR.BS_Evaluator;
import BasisSet.Atom.Gaussian.RKB.IBS_EValuator;
import BasisSet.Atom.Slater.RKB.IBS_Evaluator;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;

import qchem.BasisSet1.DB_Cache;

using json = nlohmann::json;

namespace BasisSet1::Atom
{

BasisSet1::Real_BS* Factory(Type type,const nlohmann::json& js,size_t Z)
{
    return Factory(type,js,Atom_EC(Z));
}

BasisSet1::Real_BS* Factory(const nlohmann::json& js,size_t Z)
{
    Type type=js["type"].template get<Type>();
    return Factory(type,js,Z);
}

BasisSet1::Real_BS* Factory(Type type, const nlohmann::json& js,const ElectronConfiguration& ec)
{
    if (BasisSet1::theGlobalCache==0)
        BasisSet1::theGlobalCache=new BasisSet1::IntegralsCache_RAM<double>(true);     
    size_t N=js["N"];
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    BasisSet1::Real_BS* bs=0;
    switch (type)
    {
    case Type::Slater:
    {
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet<Slater_BS>(N,emin,emax,ec);
        break;
    }
    case Type::Gaussian:
    {
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet<Gaussian_BS>(N,emin,emax,ec);
        break;
    }
    case Type::BSpline6:
    {
        double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
        bs=new BasisSet<BSpline_BS<6>>(N,rmin,rmax,ec);
        break;   
    }
//     case Type::BSpline16:
//     {
//         double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         bs=new AtomBS::BSpline1::BasisSet<6,BSpline_r_BS>(N,rmin,rmax,ec);
//         break;   
//     } 
//     case Type::BSpline9:
//     {
//         double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         bs=new AtomBS::BSpline::BasisSet<9>(N,rmin,rmax,ec);
//         break;   
//     } 
// case Type::BSpliner6:
//     {
//         double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         bs=new AtomBS::BSpline::BasisSet_r<6>(N,rmin,rmax,ec);
//         break;   
//     } 
//     case Type::BSpliner9:
//     {
//         double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         bs=new AtomBS::BSpline::BasisSet_r<9>(N,rmin,rmax,ec);
//         break;   
//     } 

    case Type::Slater_RKB:
    {
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet_RKB<Slater_IBS,Slater_RKBS_IBS>(N,emin,emax,ec);
        break;
    }
    case Type::Gaussian_RKB:
    {
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet_RKB<Gaussian_IBS,Gaussian_RKBS_IBS>(N,emin,emax,ec);
        break;
    }
//     case Type::BSpline_RKB:
//     {
//         assert(false);
//         //double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         // bs=new Atom_kappa::BSpline::BasisSet<6>(N,rmin,rmax,LMax);
//         break;   
//     } 




    }
    
    assert(bs);
    return bs;
}

} //namespace 
