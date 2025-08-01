// File: SCFParams.C Parameters used for controling SCF iteration and convergence.
export module qchem.SCFParams;
import qchem.Types;

export struct SCFParams
{
    size_t NMaxIter;         //Max allowed number of iterations
    double MinDeltaRo;       //Minimum delta in charge density for convergence.
    double MinDelE;          //Minimum delta energy convergence.
    double MinError;         //Minimum error from SCF accelerator.  i.e. [F,D] (orbital gradient).
    double StartingRelaxRo;  //relaxation for mixing Ro.  Dynamically adjusted during iterations.
    double MergeTol;         //Merge eigen levels (like Px,Py Pz) that are equal within +/- MergeTol
    bool   Verbose;          //Display iteration details.        
};


