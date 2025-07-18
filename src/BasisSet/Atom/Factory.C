// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.

#include <cassert>
#include <nlohmann/json.hpp>
#include <BasisSet/Factory.H>

#include "l/Slater_BS.H"
#include "ml/Slater_BS.H"
#include "kappa/Slater_BS.H"
#include "l/Gaussian_BS.H"
#include "ml/Gaussian_BS.H"
#include "kappa/Gaussian_BS.H"
#include "l/BSpline_BS.H"
#include "ml/BSpline_BS.H"
using json = nlohmann::json;
import qchem.Symmetry.AtomEC;

namespace BasisSetAtom
{
enum class AngularType {Yl,Ylm};

BasisSet* Factory(Type type,const nlohmann::json& js,size_t Z)
{
    return Factory(type,js,Atom_EC(Z));
}

BasisSet* Factory(const nlohmann::json& js,size_t Z)
{
    Type type=js["type"].template get<Type>();
    return Factory(type,js,Z);
}

BasisSet* Factory(Type type, const nlohmann::json& js,const ElectronConfiguration& ec)
{
    size_t N=js["N"];
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    AngularType atype = aec.IsMagnetic() ? AngularType::Ylm : AngularType::Yl;
    BasisSet* bs=0;
    switch (type)
    {
    case Type::Slater:
    {
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        switch (atype)
        {
        case AngularType::Yl:
            bs=new Atoml::Slater::BasisSet(N,emin,emax,LMax);
            break;
        
        case AngularType::Ylm:
            bs=new Atom_ml::Slater::BasisSet(N,emin,emax,ec);
        }
        break;
    }
    case Type::Gaussian:
    {
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        switch (atype)
        {
        case AngularType::Yl:
            bs=new Atoml::Gaussian::BasisSet(N,emin,emax,LMax);
            break;
        
        case AngularType::Ylm:
            bs=new Atom_ml::Gaussian::BasisSet(N,emin,emax,ec);
        }
        break;
    }
    case Type::BSpline:
    {
        double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
        switch (atype)
        {
        case AngularType::Yl:
            bs=new Atoml::BSpline::BasisSet<6>(N,rmin,rmax,LMax);
            break;
        
        case AngularType::Ylm:
            bs=new Atom_ml::BSpline::BasisSet<6>(N,rmin,rmax,ec);
        }
        break;   
    } 

    case Type::Slater_RKB:
    {
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new Atom_kappa::Slater::BasisSet(N,emin,emax,LMax);
        break;
    }
    case Type::Gaussian_RKB:
    {
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new Atom_kappa::Gaussian::BasisSet(N,emin,emax,LMax);
        break;
    }
    case Type::BSpline_RKB:
    {
        assert(false);
        //double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
        // bs=new Atom_kappa::BSpline::BasisSet<6>(N,rmin,rmax,LMax);
        break;   
    } 




    }
    
    assert(bs);
    return bs;
}

} //namespace 
