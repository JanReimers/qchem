// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.
module; 
#include <cassert>
#include <nlohmann/json.hpp>

module qchem.BasisSet1.Atom.Factory;
import qchem.Symmetry.AtomEC;
// import qchem.BasisSet.Atom.Gaussian.RKB.BS;
// import qchem.BasisSet.Atom.Slater.RKB.BS;
// import qchem.BasisSet.Atom.Gaussian.NR.BS;
// import qchem.BasisSet.Atom.Slater.NR.BS;
// import qchem.BasisSet.Atom.BSpline.NR.BS;

import qchem.BasisSet1.Atom.BSpline.NR.BS;
import BasisSet.Atom.BSpline.NR.BS_Evaluator;

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
    size_t N=js["N"];
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    BasisSet1::Real_BS* bs=0;
    switch (type)
    {
    // case Type::Slater:
    // {
    //     double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
    //     bs=new AtomBS::Slater::BasisSet(N,emin,emax,ec);
    //     break;
    // }
    // case Type::Gaussian:
    // {
    //     double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
    //     bs=new AtomBS::Gaussian::BasisSet(N,emin,emax,ec);
    //     break;
    // }
    case Type::BSpline6:
    {
        double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
        bs=new BSpline1::BasisSet<6,BSpline_r_BS>(N,rmin,rmax,ec);
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

//     case Type::Slater_RKB:
//     {
//         double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
//         bs=new AtomBS::Slater_RKB::BasisSet(N,emin,emax,LMax);
//         break;
//     }
//     case Type::Gaussian_RKB:
//     {
//         double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
//         bs=new AtomBS::Gaussian_RKB::BasisSet(N,emin,emax,LMax);
//         break;
//     }
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
