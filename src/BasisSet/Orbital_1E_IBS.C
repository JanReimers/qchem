// File: BasisSet/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron calculation.
module;

export module qchem.Orbital_1E_IBS;
export import qchem.IrrepBasisSet;
export import qchem.Cluster;

import qchem.LASolver;


//! \brief Interface for Laplacian integrals using in kinetic energy calculations.
export template <class T> class Integrals_Kinetic
{
public:
    //! Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual const SMatrix<T>& Kinetic() const=0;
};
//! \brief Interface for electron-nucleus attraction integrals.
export template <class T> class Integrals_Nuclear
{
public:
    //! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
    virtual const SMatrix<T>& Nuclear(const Cluster*) const=0;   
};

//
// Define an orbital irrep basis set which supports integrals for SCF orbital calculations.
// Mix-in the integral interfaces required for an orbital basis. 
//
export template <class T> class Orbital_IBS
    : public virtual IrrepBasisSet<T> //brings in Integrals_Overlap<T>
    , public virtual Integrals_Kinetic<T> 
    , public virtual Integrals_Nuclear<T> 
{
public:    
    virtual void         Set(const LAParams&)=0;
    virtual LASolver<T>* CreateSolver() const=0;
};

export typedef Orbital_IBS<double>    Real_OIBS;
export typedef Orbital_IBS<dcmplx> Complex_OIBS;
