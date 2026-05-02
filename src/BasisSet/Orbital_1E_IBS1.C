// File: BasisSet/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron calculation.
module;
#include <cassert>
export module qchem.Orbital_1E_IBS1;
export import qchem.IrrepBasisSet1;
export import qchem.Cluster;


//! \brief Interface for Laplacian integrals using in kinetic energy calculations.
//! Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
export template <class T> class Integrals_Kinetic1 : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T> MakeKinetic() const=0;
    
    const smat_t<T>& Kinetic() const
    {
        auto cache=theGlobalCache;
        assert(cache);
        if (!cache->Has(IntegralsCache_Base::I2C::Kinetic,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID())))
            cache->Set(MakeKinetic()); //Uses the key from the Has call.
        return cache->GetSMat(); //Uses the iterator from the Has call.
    }
};

//! \brief Interface for electron-nucleus attraction integrals.
//! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
export template <class T> class Integrals_Nuclear1 : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T> MakeNuclear(const Cluster* cl) const=0;
    
    const smat_t<T>& Nuclear(const Cluster* cl) const
    {
        assert(cl);
        auto cache=theGlobalCache;
        assert(cache);
        if (!cache->Has(IntegralsCache_Base::I2n::Nuclear,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),cl->ID()))
            cache->Set(MakeNuclear(cl)); //Uses the key from the Has call.
        return cache->GetSMat(); //Uses the iterator from the Has call.
    }
};


//
// Define an orbital irrep basis set which supports integrals for SCF orbital calculations.
// Mix-in the integral interfaces required for an orbital basis. 
//
export template <class T> class Orbital_IBS1
    : public IrrepBasisSet1<T> //brings in Integrals_Overlap<T>
    , public virtual Integrals_Kinetic1<T> 
    , public virtual Integrals_Nuclear1<T> 
{
public:    
    Orbital_IBS1(const Irrep_QNs::sym_t& sym) : IrrepBasisSet1<T>(sym) {};
};

export typedef Orbital_IBS1<double>    Real_OIBS1;
export typedef Orbital_IBS1<dcmplx> Complex_OIBS1;
