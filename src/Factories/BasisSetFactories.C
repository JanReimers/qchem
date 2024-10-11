// File: BasisSetFactories.C

#include "BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/BasisFunction.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/PolarizedGaussian/BasisFunction.H"
#include "Imp/BasisSet/PolarizedGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IntegralEngine.H"
#include "Imp/BasisSet/PolarizedGaussian/Block.H"
//#include "BasisSetImplementation/PlaneWave/PlaneWaveBF.H"
//#include "BasisSetImplementation/PlaneWave/PlaneWaveBS.H"
//#include "BasisSetImplementation/PlaneWave/PlaneWaveIE.H"

#include "DFTDataBase/HeapDB/HeapDB.H"

#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include <UnitSymmetryQN.H>
//#include "BasisSetImplementation/PlaneWave/BlochQN.H"

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
    if (Name==typeid(SphericalGaussian::BasisFunction).name()) return new SphericalGaussian::BasisFunction;
    if (Name==typeid(PolarizedGaussian::BasisFunction).name()) return new PolarizedGaussian::BasisFunction;
//    if (Name==typeid(PlaneWaveBF        ).name()) return new         PlaneWaveBF;

    std::cout << "Unknown basis function type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

namespace PolarizedGaussian
{

Block* Block::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(Block).name()) return new Block;

    std::cout << "Unknown basis function block type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

} //namespace PolarizedGaussian

//##################################################################
//
//  Basis set factory, reads in name or BasisSet derived
//  class and make a new object using the default constructor.
//

IrrepBasisSet* IrrepBasisSet::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalGaussian::IrrepBasisSet).name()) return new SphericalGaussian::IrrepBasisSet;
    if (Name==typeid(PolarizedGaussian::IrrepBasisSet).name()) return new PolarizedGaussian::IrrepBasisSet;
//    if (Name==typeid(        PlaneWaveBS).name()) return new         PlaneWaveBS;

    std::cout << "Unknown irrep basis set type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

BasisSet* BasisSet::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalGaussian::BasisSet).name()) return new SphericalGaussian::BasisSet;
    if (Name==typeid(PolarizedGaussian::BasisSet).name()) return new PolarizedGaussian::BasisSet;

    std::cout << "Unknown basis group type :" << Name << std::endl;
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

//template <class T> NumericalIE<T>* NumericalIE<T>::Factory(std::istream& is)
//{
//    std::string Name=StreamableObject::PeekAtName(is);
//    if (Name==typeid(NumericalIEImp<double>).name()) return new NumericalIEImp<double>;
//    
//    std::cout << "Unknown integral engine type :" << Name << std::endl;
//    exit(-1);
//    return NULL;
//}

//template <> NumericalIE<double>* NumericalIE<double>::Factory(std::istream& is)
//{
//  std::string Name=PeekAtName(is);
//  if (Name==typeid(NumericalIEImp<double>).name()) return new NumericalIEImp<double>;
//  
//  std::cout << "Unknown integral engine type :" << Name << std::endl;
//  exit(-1);
//  return NULL;
//}

//template <class T> AnalyticIE<T>* AnalyticIE<T>::Factory(std::istream& is)
//{
//    std::string Name=StreamableObject::PeekAtName(is);
//    if (Name==typeid(SphericalGaussianIE1).name()) return new SphericalGaussianIE1;
//    if (Name==typeid(PolarizedGaussianIE1).name()) return new PolarizedGaussianIE1;
//
//    std::cout << "Unknown integral engine type :" << Name << std::endl;
//    exit(-1);
//    return NULL;
//}
//
//template <> AnalyticIE<double>* AnalyticIE<double>::Factory(std::istream& is)
//{
//  std::string Name=PeekAtName(is);
//  if (Name==typeid(SphericalGaussianIE1).name()) return new SphericalGaussianIE1;
//  if (Name==typeid(PolarizedGaussianIE1).name()) return new PolarizedGaussianIE1;
//
//  std::cout << "Unknown integral engine type :" << Name << std::endl;
//  exit(-1);
//  return NULL;
//}

//template <> NumericalIE<std::complex<double> >* NumericalIE<std::complex<double> >::Factory(std::istream& is)
//{
//  std::string Name=PeekAtName(is);
////  if (Name==typeid(        PlaneWaveIE).name()) return new         PlaneWaveIE;
//
//  std::cout << "Unknown complex integral engine type :" << Name << std::endl;
//  exit(-1);
//  return NULL;
//}



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
//    if (Name==typeid(            BlochQN).name()) return new             BlochQN;

    std::cout << "Unknown Quantum Number type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

