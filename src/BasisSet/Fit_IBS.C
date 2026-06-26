// File: BasisSet/Fit_IBS.C  Interface for a fitting Basis Set.
module;
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.ScalarFunction;
export import qchem.Mesh1;            // qcMesh1::Mesh / MeshParams -- the fit quadrature mesh + knobs
import qchem.Structure;               // Structure (SetMesh builds the Becke mesh from it)

export namespace BasisSet
{

 //! \brief Interface for fit basis set that can all integrals required by least the squares Fitting module.
class Fit_IBS
    : public virtual Integrals_Overlap<double>
    , public virtual Real_IBS //Real Irrep basis Set
{
public:
    using Integrals_Overlap<double>::Overlap;
    using Integrals_Overlap<double>::MakeOverlap;
    const  rvec_t& Charge   () const;
    const rsmat_t& Repulsion() const;
    const  rmat_t& Repulsion(const Fit_IBS& b) const;
    const rsmat_t& InvOverlap() const;
    const rsmat_t& InvRepulsion() const;

    //! Build and OWN the fit quadrature mesh (Becke, from the structure).  Called by the
    //! CreateCD/VxcFitBasisSet creators, which already hold the Structure.
    void SetMesh(const Structure& cl, const qcMesh1::MeshParams& mp);

    // Numerical (mesh-quadrature) versions -- run over the fit basis's OWN mesh (itsMesh).
    typedef ScalarFunction<double> Sf;
    virtual const rvec_t& Norm   ()           const; //!< 1/sqrt(<f_a|f_a>), cached
    virtual rvec_t        Overlap(const Sf& f) const; //!< projection <f_a|f> (Vxc fit RHS; NOT cached)

protected:
    virtual  rvec_t MakeCharge      () const=0;
    virtual rsmat_t MakeRepulsion   () const=0;
    virtual  rmat_t MakeRepulsion   (const Fit_IBS&) const=0;
    virtual rsmat_t MakeInvOverlap  () const;
    virtual rsmat_t MakeInvRepulsion() const;

    virtual  rvec_t MakeNorm   () const; //Numerical, over itsMesh.

private:
    qcMesh1::Mesh itsMesh;   //!< the fit basis's own quadrature mesh.
};

}//namespace

