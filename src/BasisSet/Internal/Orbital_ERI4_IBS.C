// File: BasisSet/Internal/Orbital_ERI4_IBS.C  The ERI4 SUBSTRATE face behind the HF contraction.
module;
export module qchem.BasisSet.Internal.Orbital_ERI4_IBS;
export import qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Internal.ERI4;

export namespace qchem::BasisSet
{

//--------------------------------------------------------------------------------
//
//! \brief The ERI4 SUBSTRATE face: a basis that answers the HF contraction by building and caching
//! 4-index \f$(ab|cd)\f$ blocks.
//!
//! This is the implementation detail behind \c Orbital_HF_IBS, not part of it (doc/CleanupCandidates.md
//! R1.7).  It lives in an \c Internal module deliberately: NOTHING outside qcBasisSet needs an ERI4 --
//! the Hamiltonian and ChargeDensity libraries consume only the contraction face.  Its clients are the
//! concrete AO bases that implement it (Atom, Molecule, the RKB/DHF pair) plus the unit tests that pin
//! the cache and the bra-ket symmetry.
//!
//! An implementor supplies just MakeDirect()/MakeExchange(); the four \c Accumulate methods of the
//! contraction face are then implemented HERE, once, in terms of the cached Direct()/Exchange().
//! A basis with no ERI4 substrate (the SALC decorator) simply does not derive from this class.
//
template <class T> class Orbital_ERI4_IBS
    : public virtual Orbital_HF_IBS<T>
{
public:
    //! Build the Coulomb block \f$J(ab,cd)\f$.  Only called once for a given {radial,angular} ID pair.
    virtual ERI4       MakeDirect  (const Orbital_ERI4_IBS<T>& c) const=0;
    //! Build the exchange block \f$K(ab,cd)\f$.  Only called once for a given {radial,angular} ID pair.
    virtual ERI4       MakeExchange(const Orbital_ERI4_IBS<T>& c) const=0;
    //! Cached accessor for MakeDirect().  CANONICAL orientation only (min BasisSetID first) -- the
    //! cache THROWS on the partner orientation, which is what makes the bra-ket fusion below safe.
    const   ERI4&          Direct  (const Orbital_ERI4_IBS<T>& c) const;
    //! Cached accessor for MakeExchange().  Canonical orientation only, as for Direct().
    const   ERI4&          Exchange(const Orbital_ERI4_IBS<T>& c) const;

    // The contraction face, implemented over the ERI4 blocks.  See Orbital_HF_IBS for the contract;
    // the comments in Internal/Imp/Orbital_ERI4_IBS.C cover how the bra-ket fusion delivers it.
    virtual void AccumulateDirect  (rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const;
    virtual void AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const;
    virtual void AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const smat_t<T>& Di, const smat_t<T>& Dj,
                                      const Orbital_HF_IBS<T>* cd) const;
    virtual void AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const smat_t<T>& Di, const smat_t<T>& Dj,
                                        const Orbital_HF_IBS<T>* cd) const;
protected:
    //! Cross-cast a contraction-face partner down to the ERI4 substrate (abstract base to abstract
    //! base, the sanctioned direction).  THROWS, naming both bases, if the partner has no ERI4 face:
    //! an ERI4 basis and a substrate-free basis (a SALC decorator) cannot be paired in one contraction,
    //! because there is no \f$(ab|cd)\f$ block spanning the two.  Formerly this mistake produced an
    //! empty ERI4 and hence a silent zero contribution.
    const Orbital_ERI4_IBS<T>& Substrate(const Orbital_HF_IBS<T>* cd) const;
};

typedef Orbital_ERI4_IBS<double>    Real_ERI4_OIBS;
typedef Orbital_ERI4_IBS<dcmplx> Complex_ERI4_OIBS;

} //namespace
