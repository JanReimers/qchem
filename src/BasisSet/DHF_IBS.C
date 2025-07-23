// File: DHF_IBS.H  Interface for a Dirac-Hartree-Fock (HF) Orbital Irrep Basis Set.
module;

export module qchem.DHF_IBS;
export import qchem.Irrep_BS;
import qchem.BasisSet.Integrals;

export template <class T> class Orbital_RKBS_IBS;
//! \brief Interface for L-S cross kinetic matrix used in relativistic calculations.
export template <class T> class Integrals_XKinetic : public virtual Integrals_Base<T>
{
public:
    //! L/S cross Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual const Matrix<T>& Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const=0;
};

//! \brief Interface for one electron integrals used in Dirac-Hartree-Fock (DHF) calculations.
export template <class T> class Integrals_RKB
: public virtual Integrals_Overlap<T>
, public virtual Integrals_Kinetic<T>
, public virtual Integrals_Nuclear<T>
, public virtual Integrals_RestMass<T>
{
public:
    
};

//! \brief Interface for Large-Component one electron integrals used in Dirac-Hartree-Fock (DHF) calculations.
export template <class T> class Integrals_RKBL
: public virtual Integrals_Overlap<T>
, public virtual Integrals_XKinetic<T>
, public virtual Integrals_Nuclear<T>
{
public:
    
};

//! \brief Interface for Small-Component one electron integrals used in Dirac-Hartree-Fock (DHF) calculations.
export template <class T> class Integrals_RKBS 
: public virtual Integrals_Kinetic<T> //Serves as the overlap.
, public virtual Integrals_Nuclear<T>
{
public:
   
 };

 
export template <class T> class Orbital_RKB_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKB<T> 
{
public:
    
};

export template <class T> class Orbital_RKBL_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKBL<T> //One electron integrals used for everything
{
public:
    //int GetKappa() const;
};

export template <class T> class Orbital_RKBS_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_RKBS<T> //One electron integrals used for everything
{
public:
    // using IrrepBasisSet::size;
    //int GetKappa() const;
    virtual void InsertBasisFunctions(const Orbital_RKBL_IBS<T>* l)=0;
};

