// File: Hamiltonian/Imp/Factory.C  Construct and return various Hamiltonian types.
module;
#include <cassert>
#include "Hamiltonians.H"
module qchem.Hamiltonian.Factory;
namespace HamiltonianF
{
    
    Hamiltonian* Factory(Model m,Pol p, const cl_t& cl)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
            {
                switch (m)
                {
                    case Model::E1:
                        h=new Ham_1E(cl);
                        break;
                    case Model::HF:
                        h=new Ham_HF_U(cl);
                        break;
                    case Model::DE1:
                        h=new Ham_DHF_1E(cl);
                        break;
                    case Model::DHF:
                        assert(false); //DHF is always polarized?
                        h=new Ham_DHF(cl);
                        break;
                }
                break;
            }
            case Pol::Polarized:
            {
                switch (m)
                {
                case Model::E1:
                    h=new Ham_1E(cl);
                    break;
                case Model::HF:
                    h=new Ham_HF_P(cl);
                    break;
                case Model::DE1:
                    h=new Ham_DHF_1E(cl);
                    break;
                case Model::DHF:
                    h=new Ham_DHF(cl);
                    break;
                }
            break;
            }
        }
        assert(h);
        return h;
    }
    //DFT version
    Hamiltonian* Factory(Pol p,const cl_t& cl,ExFunctional* ex  , const MeshParams& mp, const BasisSet* bs)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
                h=new Ham_DFT_U(cl,ex,mp,bs);
                break;
            case Pol::Polarized:
                h=new Ham_DFT_P(cl,ex,mp,bs);
                break;
        }
        assert(h);
        return h;
    }

    Hamiltonian* Factory(Pol p,const cl_t& cl,double alpha  , const MeshParams& mp, const BasisSet* bs)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
                h=new Ham_DFT_U(cl,alpha,mp,bs);
                break;
            case Pol::Polarized:
                h=new Ham_DFT_P(cl,alpha,mp,bs);
                break;
        }
        assert(h);
        return h;
    }

}