// File: Hamiltonian/Factory.C Construct and return various Hamiltonian types.
module;
#include <memory>
export module qchem.Hamiltonian.Factory;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Types;
import qchem.Mesh;
import qchem.Structure;


export namespace qchem::Hamiltonian
{
    typedef std::shared_ptr<const Structure> st_t;
    //! The one Hamiltonian method/model enum (unified -- HF, DFT, and the relativistic/1-e cases that the
    //! GUI never uses but the unit tests rely on).  E1 = 1 electron, DE1 = Dirac 1 electron.  The DFT
    //! members name a FUNCTIONAL: Xalpha (Slater exchange, tuned alpha) and LDA (parameter-free Dirac
    //! exchange + VWN5 correlation -- the molecular default).  Next functional milestone: PBE (a GGA --
    //! needs density-gradient machinery on the mesh, not just an enum value; see doc/FacadeDFTPlan.md).
    enum class Model {E1,HF,DE1,DHF, Xalpha,LDA};
    enum class Pol   {UnPolarized,Polarized};

    //! True for the DFT members (need a mesh + fit basis + SAD seed); false for HF/1-e/Dirac.
    bool IsDFT(Model);

    Hamiltonian* Factory(Model,Pol,const st_t& st);
    Hamiltonian* Factory(Pol,const st_t& st,ExFunctional* , const qcMesh::MeshParams&, const bs_t*); //DFT version
    Hamiltonian* Factory(Pol,const st_t& st,double alpha  , const qcMesh::MeshParams&, const bs_t*); //DFT version

    //! The unified resolver: turn a Model token into the concrete polymorphic Hamiltonian in ONE switch.
    //! The DFT extras (mesh, orbital basis, xalpha) are used only for the DFT members and ignored for
    //! HF/1-e/Dirac, so a single call serves every method.  This is also the compact "default Hamiltonian"
    //! entry the unit tests want -- no manual mesh/functional/fit assembly.
    Hamiltonian* Factory(Model,Pol,const st_t& st, const qcMesh::MeshParams&, const bs_t*, double xalpha);

} // namespace
