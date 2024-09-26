// File: BasisSetFactories.C

#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBF.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianIE.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBF.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBS.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianIE.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "BasisSetImplementation/PlaneWave/PlaneWaveBF.H"
#include "BasisSetImplementation/PlaneWave/PlaneWaveBS.H"
#include "BasisSetImplementation/PlaneWave/PlaneWaveIE.H"

#include "DFTDataBase/HeapDB/HeapDB.H"

#include "BasisSetImplementation/NumericalIE.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalSymmetryQN.H"
#include "BasisSetImplementation/UnitSymmetryQN.H"
#include "BasisSetImplementation/PlaneWave/BlochQN.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>


//##################################################################
//
//  Basis function factory, reads in name or BasisFunction derived
//  class and make a new object using the default constructor.
//
BasisFunction* BasisFunction::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalGaussianBF).name()) return new SphericalGaussianBF;
    if (Name==typeid(PolarizedGaussianBF).name()) return new PolarizedGaussianBF;
    if (Name==typeid(PlaneWaveBF        ).name()) return new         PlaneWaveBF;

    std::cout << "Unknown basis function type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

BasisFunctionBlock* BasisFunctionBlock::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(BasisFunctionBlock).name()) return new BasisFunctionBlock;

    std::cout << "Unknown basis function block type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

//##################################################################
//
//  Basis set factory, reads in name or BasisSet derived
//  class and make a new object using the default constructor.
//

BasisSet* BasisSet::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalGaussianBS).name()) return new SphericalGaussianBS;
    if (Name==typeid(PolarizedGaussianBS).name()) return new PolarizedGaussianBS;
//    if (Name==typeid(        PlaneWaveBS).name()) return new         PlaneWaveBS;

    std::cout << "Unknown basis set type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

//##################################################################
//
//  Integral data base factory, reads in name of IntegralDataBase
//  derived class and makes a new object using the default constructor.
//

template <class T> IntegralDataBase<T>* IntegralDataBase<T>::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(HeapDB<T>).name()) return new HeapDB<T>;

    std::cout << "Unknown integral database type :'" << Name << "'" << std::endl;
    std::cout << "Lookgin for                    :'" << typeid(HeapDB<T>).name() << "'" << std::endl;
    exit(-1);
    return NULL;
}

template IntegralDataBase<double>* IntegralDataBase<double>::Factory(std::istream&);
template IntegralDataBase<std::complex<double> >* IntegralDataBase<std::complex<double> >::Factory(std::istream&);

/*
IntegralDataBase<double>* IntegralDataBase<double>::Factory(std::istream& is)
{
  string Name=PeekAtName(is);
  if (Name==typeid(HeapDB<double>).name()) return new HeapDB<double>;

  cout << "Unknown integral database type :'" << Name << "'" << std::endl;
	cout << "Lookgin for                    :'" << typeid(HeapDB<double>).name() << "'" << std::endl;
  exit(-1);
  return NULL;
}

IntegralDataBase<std::complex<double> >* IntegralDataBase<std::complex<double> >::Factory(std::istream& is)
{
  string Name=PeekAtName(is);
  if (Name==typeid(HeapDB<std::complex<double> >).name()) return new HeapDB<std::complex<double> >;

  cout << "Unknown integral database type :'" << Name << "'" << std::endl;
	cout << "Lookgin for                    :'" << typeid(HeapDB<std::complex<double> >).name() << "'" << std::endl;
  exit(-1);
  return NULL;
}
*/

//##################################################################
//
//  Integral engine factory, reads in name of IntegralEngine derived
//  class and makes a new object using the default constructor.
//

template <class T> IntegralEngine<T>* IntegralEngine<T>::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalGaussianIE).name()) return new SphericalGaussianIE;
    if (Name==typeid(NumericalIE<double>).name()) return new NumericalIE<double>;
    if (Name==typeid(PolarizedGaussianIE).name()) return new PolarizedGaussianIE;

    std::cout << "Unknown integral engine type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

template <> IntegralEngine<double>* IntegralEngine<double>::Factory(std::istream& is)
{
  std::string Name=PeekAtName(is);
  if (Name==typeid(SphericalGaussianIE).name()) return new SphericalGaussianIE;
  if (Name==typeid(NumericalIE<double>).name()) return new NumericalIE<double>;
  if (Name==typeid(PolarizedGaussianIE).name()) return new PolarizedGaussianIE;

  std::cout << "Unknown integral engine type :" << Name << std::endl;
  exit(-1);
  return NULL;
}

template <> IntegralEngine<std::complex<double> >* IntegralEngine<std::complex<double> >::Factory(std::istream& is)
{
  std::string Name=PeekAtName(is);
//  if (Name==typeid(        PlaneWaveIE).name()) return new         PlaneWaveIE;

  std::cout << "Unknown complex integral engine type :" << Name << std::endl;
  exit(-1);
  return NULL;
}



//##################################################################
//
//  Quantum Number factory, reads in name or BasisFunction derived
//  class and make a new object using the default constructor.
//
QuantumNumber* QuantumNumber::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalSymmetryQN).name()) return new SphericalSymmetryQN;
    if (Name==typeid(     UnitSymmetryQN).name()) return new      UnitSymmetryQN;
    if (Name==typeid(            BlochQN).name()) return new             BlochQN;

    std::cout << "Unknown Quantum Number type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

