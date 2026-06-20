// File: PolarizedGaussian1/MnD/Hermite3.C  Interface class for managing 3 function Hermite coefficients
export module qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.Hermite3;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Polarization;


export namespace BasisSet::Molecule::PolarizedGaussian1
{

//
//  Encapsulate calculation and storage of Hermite expansion coefficients
//
//   nnn
//  d
//   0
//
//  This is an interface base class for primative and contracted Hermite3 blocks.
//


class Hermite3
{
public:
    virtual ~Hermite3() {};
    virtual double operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const =0;
};

} //namespace BasisSet::Molecule::PolarizedGaussian1

