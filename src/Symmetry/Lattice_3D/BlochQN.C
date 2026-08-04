// File: Symmetry/Lattice_3D/BlochQN.C  A Quantum Number translational symmetry, i.e. a wave vector.
module;
#include <iosfwd>
export module qchem.Symmetry.Lattice_3D.BlochQN;
export import qchem.Types;
export import qchem.Symmetry;
//---------------------------------------------------------------------------------
//
//  Translational symmetry, Bloch function wave vector.
//

namespace qchem::Symmetry::Lattice_3D {
export class BlochQN : public virtual qchem::Symmetry::Symmetry
{
public:
    //! \a _N = BZ-grid divisions, \a _ik = integer grid index (0..N, used for the sequence index / identity),
    //! \a _shift = the fractional Monkhorst-Pack offset in grid steps so \f$k=(ik+shift)/N\f$: \a shift=0 is the
    //! Γ-centred grid; \a shift=½ is the classic MP offset (\f$k=\pm¼\f$ for \f$N=2\f$, i.e. CP2K's default).
    BlochQN(ivec3_t _N, ivec3_t _ik, double _weight=1.0, rvec3_t _shift={0,0,0});
    virtual size_t SequenceIndex() const;
    virtual size_t GetDegeneracy() const {return 1;} // +/-k?
    virtual size_t GetPrincipleOffset() const  {return 1;}
    virtual double GetWeight() const {return weight;} // BZ-integration weight w_k of this k-point
    //! A k-STAR of band levels is one symmetry-degenerate shell of the crystal group -- let the
    //! EnergyLevels reporting layer merge equal-eigenvalue levels across k-blocks (Symmetry doc).
    virtual bool   MergeAcrossIrreps() const {return true;}
    //! Star size of this k-point in the sampled mesh: \f$w_k N_\mathrm{mesh}\f$ (an unfolded mesh point has
    //! \f$w=1/N\f$ -> 1; an IBZ wedge representative carries its star's weight -> the member count).
    virtual size_t StarSize() const;
    virtual std::ostream&  Write(std::ostream&) const;

    rvec3_t   Getk() const {return k;}

private:
    ivec3_t N;      //This is the Brillouin zone grid size which gives context for the k vector. Used for calculating the sequence index.
    ivec3_t ik;     //Integer rep. of k.
    rvec3_t k;      //Real values.
    double  weight; //BZ-integration weight w_k (sum over the k-mesh = 1).
};

//! \brief Pry the Bloch wave vector \f$k\f$ (fractional) out of an abstract symmetry handle.
//! Throws std::bad_cast if the symmetry is not a BlochQN.  Mirrors Symmetry::Atom::Getl / Getκ:
//! an IBS constructor is handed an abstract \c sym_t and uses this helper to extract the one
//! piece of concrete information it needs (here, the crystal momentum).
export rvec3_t Getk(const sym_t&);
export rvec3_t Getk(const qchem::Symmetry::Symmetry&);
} // namespace qchem::Symmetry::Lattice_3D

