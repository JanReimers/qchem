// File: Symmetry/Angular.C Common interface for various atomic (spherical) symmetries.

export module qchem.Symmetry.Angular;
export import qchem.Symmetry;
export import qchem.Symmetry.ElectronCounts;

//---------------------------------------------------------------------------------
//
// Angular_Sym for atoms
//
export class Angular_Sym
    : public virtual Symmetry
{
public:
    virtual int GetPrincipleOffset() const {return GetL();} //Add to principle QN.  For atoms this is just l.
    // New atom specific functions.
    virtual ElCounts_l GetN(const ElCounts&) const=0;
    virtual int GetL() const=0;
};
