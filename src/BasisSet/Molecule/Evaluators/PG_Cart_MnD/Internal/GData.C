

export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.GData;
export import Common.UniqueID; 
export import qchem.Types;

export namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
{
struct GData
{
    UniqueID::IDtype  ID;
    double            Alpha; //Exponent;
    rvec3_t             R;     //Center
    int               L;     //Actually a maximum L

};

}