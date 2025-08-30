// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.
module; 
#include <cassert>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
module qchem.BasisSet.Atom.Factory;
import qchem.Symmetry.AtomEC;
import qchem.BasisSet.Atom.Gaussian.RKB.BS;
import qchem.BasisSet.Atom.Slater.RKB.BS;
import qchem.BasisSet.Atom.Gaussian.NR.BS;
import qchem.BasisSet.Atom.Slater.NR.BS;
import qchem.BasisSet.Atom.Internal.l.BSplineBS;

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
            bs=new Atoml::Slater::BasisSet(N,emin,emax,ec);
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
            bs=new Atoml::Gaussian::BasisSet(N,emin,emax,ec);
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
            bs=new Atoml::BSpline::BasisSet<6>(N,rmin,rmax,ec);
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
