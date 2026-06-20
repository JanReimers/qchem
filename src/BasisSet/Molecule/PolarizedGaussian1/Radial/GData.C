

export module qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GData;
export import Common.UniqueID; 
export import qchem.Types;

export namespace BasisSet::Molecule::PolarizedGaussian1
{
struct GData
{
    UniqueID::IDtype  ID;
    double            Alpha; //Exponent;
    rvec3_t             R;     //Center
    int               L;     //Actually a maximum L

};

}