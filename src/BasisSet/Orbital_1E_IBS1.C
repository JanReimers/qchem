// File: BasisSet/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron calculation.
module;
#include <cassert>
#include <string>
export module qchem.Orbital_1E_IBS1;
export import qchem.IrrepBasisSet1;
export import qchem.Cluster;
import qchem.BasisSet.DB_Cache1;

//  The are used for caching 1) radial Slater integrals R_k(abcd) 2) Direct/Exchange integrals
export class IrrepBasisSet_IDs
{
public:
    virtual std::string  RadialID() const=0;
    virtual std::string AngularID() const=0;
    virtual std::string Name     () const=0;

};

//--------------------------------------------------------------------------------
//
//! The method of integral evaluation is of course strongly dependant on the
//! precise details of basis functions or basis set.  
//! All integral functions except MakeNormalization return normalized integrals.
//! Interfaces for 1 electron integrals used for all Irrep basis sets: 1E,Fit,HF,DFT,DHF  
//! The calls return matrix refrences which implies they are buffered behind the scenes.
//

//! \brief Interface for overlap integrals.
//! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
export template <class T> class Integrals_Overlap1 : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T>  MakeOverlap() const=0; //Only called once for a given {radial,angular} ID pair.
    const   smat_t<T>&     Overlap() const
    {
        auto cache=theGlobalCache;
        assert(cache);
        return cache->Has(IntegralsCache_Base::I2C::Overlap,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()))
            ? cache->GetSMat() : cache->Set(MakeOverlap());
    }
};

//! \brief Interface for Kinetic energy integrals.
//! Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
export template <class T> class Integrals_Kinetic1 : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T> MakeKinetic() const=0; //Only called once for a given {radial,angular} ID pair.
    
    const smat_t<T>& Kinetic() const
    {
        auto cache=theGlobalCache;
        assert(cache);
        return cache->Has(IntegralsCache_Base::I2C::Kinetic,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()))
            ? cache->GetSMat() : cache->Set(MakeKinetic());
       
    }
};

//! \brief Interface for electron-nucleus attraction integrals.
//! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
export template <class T> class Integrals_Nuclear1 : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T> MakeNuclear(const Cluster* cl) const=0; //Only called once for a given {radial,angular, cluster} ID triple.
    
    const smat_t<T>& Nuclear(const Cluster* cl) const
    {
        assert(cl);
        auto cache=theGlobalCache;
        assert(cache);
        return cache->Has(IntegralsCache_Base::I2n::Nuclear,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),cl->ID())
            ? cache->GetSMat() : cache->Set(MakeNuclear(cl));
    }
};


//
// Define an orbital irrep basis set which supports integrals for 1-electron (1E) orbital calculations.
// Mix-in the integral interfaces required for a 1E orbital basis. 
//
export template <class T> class Orbital_1E_IBS1
    : public IrrepBasisSet1<T> //brings in symmetry and op()(r)
    , public virtual Integrals_Overlap1<T> 
    , public virtual Integrals_Kinetic1<T> 
    , public virtual Integrals_Nuclear1<T> 
{
public:    
    Orbital_1E_IBS1(const Irrep_QNs::sym_t& sym) : IrrepBasisSet1<T>(sym) {};
};

export typedef Orbital_1E_IBS1<double>    Real_OIBS1;
export typedef Orbital_1E_IBS1<dcmplx> Complex_OIBS1;
