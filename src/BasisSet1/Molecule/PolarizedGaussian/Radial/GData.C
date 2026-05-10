

export module qchem.BasisSet1.Molecule.PolarizedGaussian.Internal.GData;
export import Common.UniqueID; 
export import qchem.Types;

export namespace BasisSet1::Molecule::PolarizedGaussian
{
struct GData
{
    UniqueID::IDtype  ID;
    double            Alpha; //Exponent;
    rvec3_t             R;     //Center
    int               L;     //Actually a maximum L

};

}