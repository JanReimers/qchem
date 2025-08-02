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
    BlochQN(IVec3 _N, IVec3 _ik);
    virtual size_t SequenceIndex() const;
    virtual int GetDegeneracy() const {return 1;} // +/-k?
    virtual int GetPrincipleOffset() const  {return 1;}
    virtual std::ostream&  Write(std::ostream&) const;

    RVec3   Getk() const {return k;}

private:
    IVec3 N;  //This is the Brillouin zone grid size which gives context for the k vector. Used for calculating the sequence index.
    IVec3 ik; //Integer rep. of k.
    RVec3 k;  //Real values.
};

