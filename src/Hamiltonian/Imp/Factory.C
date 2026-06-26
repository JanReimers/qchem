// File: Hamiltonian/Imp/Factory.C  Construct and return various Hamiltonian types.
module;
#include <cassert>
module qchem.Hamiltonian.Factory;
import qchem.Hamiltonian.Internal.Hamiltonians;

namespace qchem::Hamiltonian
{
    
    Hamiltonian* Factory(Model m,Pol p, const st_t& st)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
            {
                switch (m)
                {
                    case Model::E1:
                        h=new Ham_1E(st);
                        break;
                    case Model::HF:
                        h=new Ham_HF_U(st);
                        break;
                    case Model::DE1:
                        h=new Ham_DHF_1E(st);
                        break;
                    case Model::DHF:
                        h=new Ham_DHF_U(st);
                        break;
                }
                break;
            }
            case Pol::Polarized:
            {
                switch (m)
                {
                case Model::E1:
                    h=new Ham_1E(st);
                    break;
                case Model::HF:
                    h=new Ham_HF_P(st);
                    break;
                case Model::DE1:
                    h=new Ham_DHF_1E(st);
                    break;
                case Model::DHF:
                    h=new Ham_DHF_P(st);
                    break;
                }
            break;
            }
        }
        assert(h);
        return h;
    }
    //DFT version
    Hamiltonian* Factory(Pol p,const st_t& st,ExFunctional* ex  , const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
                h=new Ham_DFT_U(st,ex,mp,bs);
                break;
            case Pol::Polarized:
                h=new Ham_DFT_P(st,ex,mp,bs);
                break;
        }
        assert(h);
        return h;
    }

    Hamiltonian* Factory(Pol p,const st_t& st,double alpha  , const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
                h=new Ham_DFT_U(st,alpha,mp,bs);
                break;
            case Pol::Polarized:
                h=new Ham_DFT_P(st,alpha,mp,bs);
                break;
        }
        assert(h);
        return h;
    }

}