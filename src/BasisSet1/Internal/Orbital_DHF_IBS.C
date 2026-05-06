// File: BasisSet1/Orbital_DHF_IBS.C Interface for a Dirac-Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
export module qchem.BasisSet1.Internal.Orbital_DHF_IBS;
export import qchem.BasisSet1.IrrepBasisSet;
export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.Cluster;

export namespace BasisSet1
{

template <class T> class Orbital_RKBS_IBS;
//  Large portion of an RKB irrep basis set for DHF calculations.
template <class T> class Orbital_RKBL_IBS
    : public virtual IrrepBasisSet<T>
    , public virtual Integrals_Overlap<T>
    , public virtual Integrals_Nuclear<T>
    , public virtual IrrepBasisSet_IDs //avoid using statements for RadialID,AngularID
{
public:
    //! L/S cross Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual const mat_t<T>&     Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    virtual       mat_t<T>  MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const=0;

};

//  Small portion of an RKB irrep basis set for DHF calculations.
template <class T> class Orbital_RKBS_IBS
    : public virtual IrrepBasisSet<T>
    , public virtual Integrals_Kinetic<T> //Serves as the overlap.
    , public virtual Integrals_Nuclear<T>
    , public virtual IrrepBasisSet_IDs //avoid using statements for RadialID,AngularID
{
public:
    virtual void Insert(const Orbital_RKBL_IBS<T>* l)=0;
};

} //namespace