

export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GData;
export import Common.UniqueID; 
export import qchem.Types;

export namespace PolarizedGaussian
{
struct GData
{
    UniqueID::IDtype  ID;
    double            Alpha; //Exponent;
    RVec3             R;     //Center
    int               L;     //Actually a maximum L

};

}