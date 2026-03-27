// File: TOrbitals.C  
module;
#include <iosfwd>
export module qchem.Orbitals.Internal.OrbitalsImp;
export import qchem.Orbitals;
export import qchem.IrrepBasisSet;
export import qchem.Symmetry.Irrep;
export import qchem.Types;



export template <class T> class TOrbitalsImp
    : public virtual Orbitals
    , public virtual TOrbitals<T>
{
    typedef typename TOrbitals<T>::ds_t ds_t; //{double,smat_t}
public:
    TOrbitalsImp(const Orbital_IBS<T>*, Spin s);
    virtual ~TOrbitalsImp();


    virtual ds_t      TakeElectrons      (double ne      )      ;
    virtual size_t    GetNumOrbitals     (               ) const;
    virtual size_t    GetNumOccOrbitals  (               ) const;
    virtual double    GetEigenValueChange(const Orbitals&) const;
    virtual DM_CD*    GetChargeDensity   () const;
    virtual void      UpdateOrbitals     (const mat_t<T>& U, const mat_t<T>& UPrime, const rvec_t& e);
    virtual Irrep_QNs GetQNs() const;
    
    virtual vec_t    <T> operator() (const RVec3&) const;
    virtual vec3vec_t<T> Gradient   (const RVec3&) const;


    virtual const_iterator begin() const {return itsOrbitals.begin();}
    virtual const_iterator end  () const {return itsOrbitals.end  ();} 
    virtual       iterator begin()       {return itsOrbitals.begin();}
    virtual       iterator end  ()       {return itsOrbitals.end  ();} 


    virtual std::ostream&          Write(std::ostream&) const;

private:
    TOrbitalsImp(const TOrbitalsImp&);

    const Orbital_IBS<T>*  itsBasisSet;
    ov_t                   itsOrbitals;
    Irrep_QNs              itsQNs;
    smat_t<T>              itsD; // DPrime=C'*Cd',  U*D*Ud, D=C*Cd (outer product)
};

