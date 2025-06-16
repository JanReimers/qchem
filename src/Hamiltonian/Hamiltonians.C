// File:: Hamiltonians.S  Create fully implemented Hamiltonians

#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Hamiltonian/Vee.H"
#include "Imp/Hamiltonian/Vxc.H"
#include <BasisSet/BasisSet.H>
#include <BasisSet/Fit_IBS.H>
#include <memory>

Ham_1E::Ham_1E(const cl_t& cl) 
{
    InsertStandardTerms(cl);
}

Ham_HF_U::Ham_HF_U(const cl_t& cl) 
{
    InsertStandardTerms(cl);
    Add(new Vee);
    Add(new Vxc);
}

#include "Imp/Hamiltonian/SlaterExchange.H"
#include "Imp/Hamiltonian/FittedVxc.H"
#include "Imp/Hamiltonian/ExchangeFunctional.H" 
#include <Cluster.H>
Ham_SHF_U::Ham_SHF_U(const cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_SHF_U(cl,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_SHF_U::Ham_SHF_U(const cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    
    FittedVxc::ex_t   XcFunct(ex);
    FittedVxc::bs_t   XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    FittedVxc::mesh_t m(cl->CreateMesh(mp));
    Add(new FittedVxc(XFitBasis, XcFunct,m));
}

#include "Imp/Hamiltonian/FittedVee.H"
#include <Cluster.H>

Ham_DFT_U::Ham_DFT_U(const cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_DFT_U(cl,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_DFT_U::Ham_DFT_U(const cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
       
    FittedVee::bs_t   CFitBasis(bs->CreateCDFitBasisSet(cl.get()));
    FittedVee::mesh_t m(cl->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,cl->GetNumElectrons()));

    FittedVxc::ex_t XcFunct(ex);
    FittedVxc::bs_t XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    Add(new FittedVxc(XFitBasis, XcFunct,m));
}

#include "Imp/Hamiltonian/VxcPol.H"
Ham_HF_P::Ham_HF_P(const cl_t& cl)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    Add(new VxcPol);
}

#include "Imp/Hamiltonian/FittedVxcPol.H"

Ham_SHF_P::Ham_SHF_P(const cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_SHF_P(cl,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_SHF_P::Ham_SHF_P(const cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
    Add(new Vee);
    
    FittedVxcPol::ex_t XcFunct(ex);
    FittedVxcPol::bs_t XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    FittedVxcPol::mesh_t  m(cl->CreateMesh(mp));
    Add(new FittedVxcPol(XFitBasis, XcFunct,m));
}

Ham_DFT_P::Ham_DFT_P(const cl_t& cl,double alpha_ex, const MeshParams& mp, const BasisSet* bs)
    : Ham_DFT_P(cl,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_DFT_P::Ham_DFT_P(const cl_t& cl,ExFunctional* ex, const MeshParams& mp, const BasisSet* bs)
{
    InsertStandardTerms(cl);
    FittedVee::bs_t CFitBasis(bs->CreateCDFitBasisSet(cl.get()));
    FittedVee::mesh_t  m(cl->CreateMesh(mp));
    Add(new FittedVee(CFitBasis,m,cl->GetNumElectrons()));
    
    FittedVxcPol::ex_t XcFunct(ex);
    FittedVxcPol::bs_t XFitBasis(bs->CreateVxcFitBasisSet(cl.get()));
    Add(new FittedVxcPol(XFitBasis, XcFunct,m));
    
}

#include "Imp/Hamiltonian/DiracKinetic.H"
#include "Imp/Hamiltonian/RestMass.H"
#include "Imp/Hamiltonian/Vnn.H"
#include "Imp/Hamiltonian/Ven.H"

Ham_DHF_1E::Ham_DHF_1E(const cl_t& cl)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(cl));
}

Ham_DHF::Ham_DHF(const cl_t& cl)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(cl));
    Add(new Ven(cl));
    // Add(new Vee());
    // Add(new VxcPol());
}