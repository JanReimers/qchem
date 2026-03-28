// File: Symmetry/BlochQN.C  A Quantum Number translational symmetry, i.e. a wave vector.
module;
#include <iosfwd>
export module qchem.Symmetry.BlochQN;
export import qchem.Types;
export import qchem.Symmetry;
//---------------------------------------------------------------------------------
//
//  Translational symmetry, Bloch function wave vector.
//

export class BlochQN : public virtual Symmetry
{
public:
    BlochQN(ivec3_t _N, ivec3_t _ik);
    virtual size_t SequenceIndex() const;
    virtual int GetDegeneracy() const {return 1;} // +/-k?
    virtual int GetPrincipleOffset() const  {return 1;}
    virtual std::ostream&  Write(std::ostream&) const;

    rvec3_t   Getk() const {return k;}

private:
    ivec3_t N;  //This is the Brillouin zone grid size which gives context for the k vector. Used for calculating the sequence index.
    ivec3_t ik; //Integer rep. of k.
    rvec3_t k;  //Real values.
};

