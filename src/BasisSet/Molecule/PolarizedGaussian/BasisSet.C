// File BasisSet/Molecule/PolarizedGaussian/BasisSet.C
module;
#include <vector>
#include <memory>
export module qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian.Reader;

import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Cluster;
import qchem.Types;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Orbital_HF_IBS;

import qchem.BasisSet.Internal.IntegralEnums;

export namespace BasisSet::Molecule::PolarizedGaussian
{

rsmat_t MakeIntegrals(IType,const PGData* ab,const Cluster*cl =0);

class IrrepBasisSet
        : public virtual Real_IBS,
          public IrrepBasisSetImp<double>,
          public PGData
    {
    public:
        typedef std::vector<std::unique_ptr<Block>> bv_t;
    
        IrrepBasisSet(Reader *, const Cluster *);
        IrrepBasisSet(const rvec_t &exponents, size_t L, const Cluster *);
        IrrepBasisSet(const rvec_t &exponents, size_t L);
        virtual ~IrrepBasisSet(); //g++ 15.2 BUG Compiler generated, or inline destructor does instance std::vector templates destructor.

        virtual size_t  GetNumFunctions() const {return PGData::size();}
        virtual size_t  size() const {return PGData::size();}
        virtual rvec_t     operator() (const rvec3_t&) const;
        virtual rvec3vec_t Gradient   (const rvec3_t&) const;

        virtual std::string RadialID () const {return PGData::RadialID();}
        virtual std::string AngularID() const {return PGData::AngularID();}
        virtual std::string BasisSetID() const {return PGData::BasisSetID();} // geometry-aware (override radial|angular default)
        virtual std::string Name     () const;
        virtual std::ostream &Write(std::ostream &) const;

    private:
        
        IrrepBasisSet(const IrrepBasisSet *bs, const bv_t&);
    
        bv_t itsBlocks;
    protected:
    };
class Orbital_IBS
    : public virtual Orbital_1E_IBS<double>
    , public virtual Orbital_HF_IBS<double>
    , public virtual Orbital_DFT_IBS<double>
    , public IrrepBasisSet
{
public:
    Orbital_IBS(Reader *, const Cluster *);
    Orbital_IBS(const rvec_t& exponents, size_t L, const Cluster *);
    Orbital_IBS(const rvec_t& exponents, size_t L);

    virtual Fit_IBS* CreateCDFitBasisSet(const Cluster *) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster *) const;

    virtual rsmat_t      MakeOverlap() const {return MakeIntegrals(PolarizedGaussian::Overlap2C,this);}
    virtual rsmat_t      MakeKinetic() const {return MakeIntegrals(Grad2,this);}
    virtual rsmat_t      MakeNuclear(const Cluster* cl) const {return MakeIntegrals(PolarizedGaussian::Nuclear,this,cl);}
    virtual ERI3<double> MakeOverlap3C  (const Fit_IBS& c) const; //Used for DFT
    virtual ERI3<double> MakeRepulsion3C(const Fit_IBS& c) const; //Used for DFT
    virtual ERI4         MakeDirect     (const ::BasisSet::Orbital_HF_IBS<double>& c) const;
    virtual ERI4         MakeExchange   (const ::BasisSet::Orbital_HF_IBS<double>& b) const;
private:
    rsmat_t Integrate(qchem::IType3C type , const GaussianRF* rc, const Polarization& pc) const;
};
// Use E prefix to avoid name clash with the interface class Fit_IBS
class EFit_IBS
    : public virtual Fit_IBS 
    , public IrrepBasisSet
{
public:
    EFit_IBS(Reader *, const Cluster *);

    virtual rsmat_t MakeOverlap() const {return MakeIntegrals(PolarizedGaussian::Overlap2C,this);}
    virtual  rvec_t MakeCharge   () const;
    virtual rsmat_t MakeRepulsion() const {return MakeIntegrals(Repulsion2C,this);}
    virtual  rmat_t MakeRepulsion(const Fit_IBS& b) const;
};
class BasisSet 
    : public virtual ::BasisSet::BasisSet<double>
    , public ::BasisSet::BasisSetImp<double>
{
public:
    BasisSet() {};
    BasisSet( Reader*, const Cluster*);
    virtual void Insert(bs_t* bs);

};

} //namespace


