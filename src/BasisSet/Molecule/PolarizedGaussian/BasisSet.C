// File PolarizedGaussian/BasisSet.C
module;
#include <vector>
#include <memory>
export module qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IEClient;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IntegralEngine;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;
import qchem.BasisSet.Molecule.PolarizedGaussian.Reader;

import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.IEClient;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.IBS_Common;
import qchem.Cluster;
import qchem.Types;
import qchem.DFT_IBS;
import qchem.HF_IBS;

export namespace PolarizedGaussian
{

class BasisFunction : public Real_BF
{
public:
    BasisFunction(                                          );
    BasisFunction(const RadialFunction*, const Polarization&, double norm);

    virtual bool   operator==(const ::Real_BF&) const;
 
    virtual std::ostream&       Write(std::ostream&   ) const;
    virtual std::istream&       Read (std::istream&   )      ;
    virtual BasisFunction* Clone(           ) const;

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;

private:
    void Insert  (const RadialFunction* theRF,const Polarization& thePol);

    const RadialFunction*  itsRadial;
    Polarization           itsPol;
    double                 itsNormalization;
};
class IrrepBasisSet
        : public virtual ::IrrepBasisSet,
          public IBS_Common,
          public IrrepIEClient
    {
    public:
        typedef std::vector<std::unique_ptr<Block>> bv_t;
    
        IrrepBasisSet(Reader *, const Cluster *);
        IrrepBasisSet(const RVec &exponents, size_t L, const Cluster *);
        IrrepBasisSet(const RVec &exponents, size_t L);

        virtual std::ostream &Write(std::ostream &) const;

    private:
        
        IrrepBasisSet(const IrrepBasisSet *bs, const bv_t&);
        void MakeBasisFunctions(const RVec &norms);

        bv_t itsBlocks;
    };
class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>,
        public virtual TOrbital_DFT_IBS<double>,
        public IrrepBasisSet,
        public Orbital_IBS_Common<double>,
        public Orbital_DFT_IBS_Common<double>,
        public Orbital_HF_IBS_Common<double>,
        public Orbital_IE

{
    typedef DB_BS_2E<double> db_t;
public:
    Orbital_IBS(const db_t* db, Reader *, const Cluster *);
    Orbital_IBS(const db_t* db, const Vector<double>& exponents, size_t L, const Cluster *);
    Orbital_IBS(const db_t* db, const Vector<double>& exponents, size_t L);

    virtual ::Fit_IBS *CreateCDFitBasisSet(const ::BasisSet*,const Cluster *) const;
    virtual ::Fit_IBS *CreateVxcFitBasisSet(const ::BasisSet*,const Cluster *) const;
    virtual IrrepBasisSet *Clone(const RVec3 &) const;
};
class Fit_IBS
    : public virtual ::Fit_IBS
    , public virtual FitIntegrals
    , public IrrepBasisSet
    , public TIBS_Common<double>
    , public Fit_IBS_Common
    , public Fit_IE
{
public:
    Fit_IBS(const DB_cache<double>* db , Reader *, const Cluster *);

    virtual ::Fit_IBS *Clone(const RVec3 &) const;
};
class BasisSet 
    : public BS_Common
    , public DB_BS_2E<double>
{
public:
    BasisSet() {};
    BasisSet( Reader*, const Cluster*);
    virtual void Insert(bs_t* bs);


    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const;
private:
    mutable CDCache cache; //Cache of all Gaussian pair charge distributions.
};

} //namespace


