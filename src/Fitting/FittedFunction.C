// File: FittedFunction.C  Linear fitted function interface
module;

export module qchem.FittedFunction;
import qchem.ScalarFunction;
import qchem.FittedFunctionClient;

//----------------------------------------------------------------------------------
//
//!  Interface for a real space function that can be fit with a basis set.\n
//! \f$ F\left(\vec{r}\right)\sim\sum_{a}c_{a}f_{a}\left(\vec{r}\right) \f$ \n
//! The expansion coefficients are calculated by least squares:\n
//! \f$ \frac{\partial}{\partial c_{b}}\int d^{3}\vec{r}\left[F\left(\vec{r}\right)-\sum_{a}c_{a}f_{a}\left(\vec{r}\right)\right]^{2}=0 \f$ \n
//! \f$ \int d^{3}\vec{r}F\left(\vec{r}\right)f_{b}\left(\vec{r}\right)=F_{b}=\sum_{a}c_{a}\int d^{3}\vec{r}f_{a}\left(\vec{r}\right)f_{b}\left(\vec{r}\right) \f$ 
//! or \f$ F_{b}=\sum_{a}c_{a}S_{ab} \f$ \n
//!  The DoFit member function calculates the expansion coefficients as: \f$ c_{b}=\sum_{a}F_{a}S_{ab}^{-1} \f$ 
//!  where \f$ S_{ab}^{-1} \f$ is the Penrose inverse of the fit basis overlap matrix. The Pensrose inverse avoids problems
//!  with linear dependensices in the fit basis.  
//!  
//!  The GetFunctionOverlap is supplied by the client (derived) class to calculate a vector of 
//! \f$ \int d^{3}\vec{r}F\left(\vec{r}\right)f_{b}\left(\vec{r}\right)=F_{b} \f$ values.
//! 
export class FittedFunction
    : public virtual ScalarFunction<double>
{
public:
    //! Find \f$ c_{b}=\sum_{a}F_{a}S_{ab}^{-1} \f$ , where \f$ F\left(\vec{r}\right) \f$ is evaluated numerically.
    virtual double DoFit(const ScalarFFClient&)=0;
    //! For the case where \f$ F\left(\vec{r}\right)=\sum_{ab}b_{a}\left(\vec{r}\right)b_{b}\left(\vec{r}\right)D_{ab} \f$ is calculated from a density matrix.
    virtual double DoFit(const DensityFFClient&)=0;
    //! do \f$ c_{b}*=factor \f$
    virtual void   ReScale         (double factor )=0; 
    //! do \f$ c_{a}=(1-f)c_{a} + fg_{a} \f$.  Assumes same basis set, and just mixes the coefficients.
    virtual void   FitMixIn        (const FittedFunction& g,double f)=0; 
    //! returns \f$ \sum_{a}\left(c_{a}-g_{a}\right)\int d^{3}\vec{r}f_{a}\left(\vec{r}\right) \f$
    virtual double FitGetChangeFrom(const FittedFunction& g) const=0;
};

