// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.

#include "Symmetry/Atom_EC.H"
#include <Factory.H>
#include "Imp/BasisSet/Atom/l/Slater_BS.H"
#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/ml/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/l/BSpline_BS.H"
#include "Imp/BasisSet/Atom/ml/BSpline_BS.H"
#include <cassert>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace BasisSetAtom
{
    enum class AngularType {Yl,Ylm};

    BasisSet* Factory(Type type,const nlohmann::json& js,size_t Z)
    {
        return Factory(type,js,Atom_EC(Z));
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
        }
    
    assert(bs);
    return bs;
    }
}
