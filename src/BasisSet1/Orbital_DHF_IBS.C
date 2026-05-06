// File: BasisSet1/Orbital_DHF_IBS.C Interface for a Dirac-Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
export module qchem.BasisSet1.Orbital_DHF_IBS;
export import qchem.BasisSet1.IrrepBasisSet;
export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.Cluster;

export namespace BasisSet1
{

//
//  Restricted Kinetic Balance (RKB) irrep basis set for DHF calculations.
//  Under the hood this will be implemented using large and small blocks.
//
template <class T> class Orbital_RKB_IBS
    : public virtual Orbital_1E_IBS<T>
    , public virtual IrrepBasisSet_IDs //avoid using statements for RadialID,AngularID
{
public:
    //! Rest mass \f$ \left\langle a\left|\left(\beta-\alpha\right)c^{2}\right|b\right\rangle =\left(\beta-\alpha\right)c^{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$
    virtual const smat_t<T>&     RestMass() const;   
    virtual       smat_t<T>  MakeRestMass() const=0;   

};


} //namespace