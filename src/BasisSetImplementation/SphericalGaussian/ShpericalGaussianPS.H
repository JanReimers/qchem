#ifndef SPERICALGAUUSIANPS_H
#define SPERICALGAUUSIANPS_H

//
//  Parameters for a spherical Gaussian basis set.
//
struct ShpericalGaussianPS
{
    SpericalGauusianPS(size_t N) :  alphas(N), Ls(N), norms(N) {};
    Vector<double> alphas; //Exponents
    Vector<size_t> Ls;     //Angular momentums
    Vector<double> norms;  //Normalization constants.
};

//class IEClientSG : public IECLient<ShpericalGaussianPS>
//{
//    
//};

#endif // SPERICALGAUUSIANPS_H
