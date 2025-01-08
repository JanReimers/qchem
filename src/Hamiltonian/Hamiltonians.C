// File:: Hamiltonians.S  Create fully implemented Hamiltonians

#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Hamiltonian/Vee.H"
#include "Imp/Hamiltonian/Vxc.H"
#include <memory>

Ham_HF_U::Ham_HF_U(cl_t& cl) 
{
    InsertStandardTerms(cl);
    Add(new Vee);
    Add(new Vxc);
}

#include "Imp/Hamiltonian/SlaterExchange.H"
#include "Imp/Hamiltonian/FittedVxc.H"
#include "Imp/Hamiltonian/ExchangeFunctional.H" 
#include <Cluster.H>
Ham_SHF_U::Ham_SHF_U(cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_SHF_U(cl,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_SHF_U::Ham_SHF_U(cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    
    std::shared_ptr<ExFunctional> XcFunct(ex);
    std::shared_ptr<const IrrepBasisSet> XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    std::shared_ptr<const Mesh>  m(cl->CreateMesh(mp));
    Add(new FittedVxc(XFitBasis, XcFunct,m));
}

#include "Imp/Hamiltonian/FittedVee.H"
#include <Cluster.H>

Ham_DFT_U::Ham_DFT_U(cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_DFT_U(cl,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_DFT_U::Ham_DFT_U(cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
       
    std::shared_ptr<const IrrepBasisSet> CFitBasis(bs->CreateCDFitBasisSet(cl.get()));
    std::shared_ptr<const Mesh>          m(cl->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,cl->GetNumElectrons()));

    std::shared_ptr<ExFunctional>        XcFunct(ex);
    std::shared_ptr<const IrrepBasisSet> XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    Add(new FittedVxc(XFitBasis, XcFunct,m));
}

#include "Imp/Hamiltonian/VxcPol.H"
Ham_HF_P::Ham_HF_P(cl_t& cl)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    Add(new VxcPol);
}

#include "Imp/Hamiltonian/FittedVxcPol.H"

Ham_SHF_P::Ham_SHF_P(cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_SHF_P(cl,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_SHF_P::Ham_SHF_P(cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    
    std::shared_ptr<ExFunctional> XcFunct(ex);
    std::shared_ptr<const IrrepBasisSet> XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    std::shared_ptr<const Mesh>  m(cl->CreateMesh(mp));
    Add(new FittedVxcPol(XFitBasis, XcFunct,m));
}

Ham_DFT_P::Ham_DFT_P(cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_DFT_P(cl,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_DFT_P::Ham_DFT_P(cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
    std::shared_ptr<const IrrepBasisSet> CFitBasis(bs->CreateCDFitBasisSet(cl.get()));
    std::shared_ptr<const Mesh>  m(cl->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,cl->GetNumElectrons()));
    
    std::shared_ptr<ExFunctional> XcFunct(ex);
    std::shared_ptr<const IrrepBasisSet> XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    Add(new FittedVxcPol(XFitBasis, XcFunct,m));
    
}

#include "Imp/Hamiltonian/DiracKinetic.H"
#include "Imp/Hamiltonian/RestMass.H"
#include "Imp/Hamiltonian/Vnn.H"
#include "Imp/Hamiltonian/Ven.H"

Ham_DHF::Ham_DHF(cl_t& cl)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    Add(new Vnn(cl));
    Add(new Ven(cl));
    // Add(new DiracVee);
    // Add(new DiracVxc);
}