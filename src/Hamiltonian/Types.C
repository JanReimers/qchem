// File: Hamiltonian/Types.C  Define types used throughout the Hamiltonian module.
module;
#include <vector>

export module qchem.Hamiltonian.Types;

export import qchem.BasisSet;
export import qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Orbital_DHF_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
import qchem.Types;   // dcmplx (for cobs_t)


export namespace qchem::Hamiltonian
{
    using bs_t    =BasisSet::BasisSet<double>;
    using fbs_t   =BasisSet::Fit_IBS;
    template <class T> using tobs_t=BasisSet::Orbital_1E_IBS<T>;  // T-parametric orbital basis
    // r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t); bare obs_t transitional (= robs_t).
    using robs_t  =tobs_t<double>;  using cobs_t=tobs_t<dcmplx>;
    using obs_t   =robs_t;
    using ohfbs_t =BasisSet::Orbital_HF_IBS<double>;
    using odftbs_t=BasisSet::Orbital_DFT_IBS<double>;
    using orkbbs_t=BasisSet::Orbital_RKB_IBS<double>;

    //! Cross-irrep view of the whole SCF Fock build, assembled by CompositeWF once per iteration and
    //! threaded through the dynamic (rho-dependent) Hamiltonian terms.  Static terms never receive it,
    //! and a dynamic term is free to ignore it (the default \c GetMatrix overload does exactly that).
    //! Its reason to exist: a term CAN exploit structure that spans more than its own irrep -- e.g. the
    //! Coulomb/exchange term banking the ERI4 bra-ket symmetry \f$J(i,j)=J(j,i)^\mathsf{T}\f$ across
    //! canonical irrep pairs (see doc/ERI4Rework.md \S5.4).  A plain struct for now; a polymorphic
    //! "answers abstract questions" type may replace it once a second consumer appears.
    template <class T> struct tHamiltonianContext
    {
        std::vector<const tobs_t<T>*> irrepBases;   //!< every participating irrep's orbital basis
    };
    using rHamiltonianContext=tHamiltonianContext<double>;  using cHamiltonianContext=tHamiltonianContext<dcmplx>;
    using HamiltonianContext =rHamiltonianContext;
}
