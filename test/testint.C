// File: TestInt.C  Test integral routines.

#include "DFTDataBase/DFTDataBase.H" 
#include "BasisSet/TBasisSet.H"
#include "BasisSet/TBasisFunction.H"
#include "BasisSet/TBasisSetBrowser.H"
#include "Cluster/Cluster.H"
#include "DFTDataBase/DFTDataBase.H"
#include "Misc/MatrixList.H"

#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/list.h"
#include "oml/stopw.h"

#include <iostream.h>
#include <fstream.h>
#include <string>
#include <GetOpt.h>
#include <math.h>
#include "Misc/Unpickle.H"
//#include "Misc/MemoryDebug.H"

template bool UnPickle(BasisSet *&,const char* ,const char*); 
template bool UnPickle(Cluster *& ,const char* ,const char*);

extern template class Vector<double>;
extern template class Matrix<double>;
extern template class SMatrix<double>;

double RelError(const Vector<double>&, const Vector<double>&);
double RelError(const Matrix<double>&, const Matrix<double>&);
double RelError(const SMatrix<double>&, const SMatrix<double>&);

//############################################################################
//
//
int main(int argc, char **argv)
{
//  MemDB m; 
  {
    

  string usage="-b<basis set> -n<numerical basis set> -c<cluster> -D<database> -v<verbose mode> -r<do repulsions>";
  string line="-------------------------------------------------------------------------";
  if(argc==1) 
  {
    cerr << "usage: " << argv[0] << usage  << endl;
    return -1;
  }

  string dfile,bfile,nfile,cfile;

  BasisSet *bs1, *bsnum1;
  TBasisSet<double> *bs, *bsnum;
  Cluster     *cl=0;
  IntegralDataBase<double>  *db ,*dbnum;
  bool verbose=false, repulsions=false;
    
  Vector<double>    an,num;
  Matrix<double>    anm,numm;
  SMatrix<double>    ansm,numsm;

  GetOpt opt (argc, argv, "b:n:c:D:vr");
  int opchar;
  
  while ((opchar = opt ()) != EOF)
    switch (opchar)
    {
      case 'b': bfile = opt.optarg; break;
      case 'n': nfile = opt.optarg; break;
      case 'c': cfile = opt.optarg; break;
      case 'D': dfile = opt.optarg; break;
      case 'v': verbose   = true  ; break;
      case 'r': repulsions= true  ; break;
      case '?': cerr << "usage: " << argv[0] << usage << endl;
                return -1;
    }

  DFTDataBase theDataBase(dfile);
  
  if (!UnPickle(bs1   ,bfile,"basis set"          ))  goto cleanup;
  if (!UnPickle(bsnum1,nfile,"numerical basis set"))  goto cleanup;
  if (!UnPickle(cl    ,cfile,"cluster"            ))  goto cleanup;
  
//	bs   =dynamic_cast<TBasisSet<double>*>(bs1   );
//	bsnum=dynamic_cast<TBasisSet<double>*>(bsnum1);
 bs = bs1->Cast(double(0));
 bsnum = bsnum1->Cast(double(0));
	assert(bs);
	assert(bsnum1);
  db   =bs   ->GetDataBase();
  dbnum=bsnum->GetDataBase();

  StreamableObject::SetOutputMode(StreamableObject::pretty);  
  cout.setf(ios::fixed,ios::floatfield);  
  cout.precision(6);  
  cout.width(10);  

//---------------------------------------------------------------------------------------------------------------
//
//  1 center integrals
//
  cout << endl << line << endl << endl;
  {  
    an =db   ->GetNormalization();
    StopWatch S;S.Start();
    num=dbnum->GetNormalization();
    S.Stop();
    if (verbose) cout << "Analytic " << an << endl << "Numerical " << num << endl;
    cout << "Normalization      RMS error : " << RelError(an,num) << " time: " << S.GetTime() << "(s)" << endl;
  }
    
  cout << endl << line << endl << endl;
  {  
    an =db   ->GetCharge();
    StopWatch S;S.Start(); 
    num=dbnum->GetCharge();
    S.Stop();
    if (verbose) cout << "Analytic " << an << endl << "Numerical " << num << endl;
    cout << "Charge             RMS error : " << RelError(an,num) << " time: " << S.GetTime() << "(s)" << endl;
  }
//---------------------------------------------------------------------------------------------------------------
//
//  2 and 3 center overlap integrals
//

    cout << endl << line << endl << endl;
  {  
    ansm =db   ->GetOverlap();
    StopWatch S;S.Start();
    numsm=dbnum->GetOverlap();
    S.Stop();
    if (verbose) cout << "Analytic " << ansm << endl << "Numerical " << numsm << endl;
    cout << "Overlap            RMS error : " << RelError(ansm,numsm) << " time: " << S.GetTime() << "(s)" << endl;
    anm=ansm;
    cout << "Eigen values of the overlap " << anm.Diagonalize();
  }
 
  cout << endl << line << endl << endl;
  {  
    anm =db   ->GetOverlap();
    StopWatch S;S.Start();
    numm=dbnum->GetOverlap(*bs);
    S.Stop();
    if (verbose) cout << "Analytic " << anm << endl << "Numerical " << numm << endl;
    cout << "Cross Overlap RMS error : " << RelError(anm,numm) << " time: " << S.GetTime() << "(s)" << endl;
  }

  cout << endl << line << endl << endl;
  {
    TBasisSetBrowser<double> bsbr(*bs);
    const ScalarFunction<double>& rsf(*bsbr);
    an =db   ->GetOverlap(rsf);
    StopWatch S;S.Start();
    num=dbnum->GetOverlap(rsf);
    S.Stop();
    if (verbose) cout << "Analytic " << an << endl << "Numerical " << num << endl;
    cout << "Overlap with scalar Function RMS error : " << RelError(an,num) << " time: " << S.GetTime() << "(s)" << endl;
  }
  cout << endl << line << endl << endl;
  {
    StopWatch S;S.Start();
    const MatrixList<double>& mlist  =db   ->GetOverlap3C(*bs);
    S.Stop();
    cout << "Overlap with Basis Set time: " << S.GetTime() << "(s)" << endl;
    S.Start();
    const MatrixList<double>& mlistnum  =dbnum   ->GetOverlap3C(*bs);
    S.Stop();
    cout << "Numerial Overlap with Basis Set time: " << S.GetTime() << "(s)" << endl;
    
    for (int i=0;i<bs->GetNumFunctions();i++)
    {
			if(verbose) cout << mlist[i] << mlistnum[i] << endl;
      cout << "Overlap with Basis Set RMS error : " << RelError(mlist[i],mlistnum[i]) << endl;
    }
  }
  cout << endl << line << endl << endl;
  
/*  cout << endl << line << endl << endl;
  {
    StopWatch S;S.Start();
    MatrixList anm_list  =db   ->GetOverlap3C(*bs);
    S.Stop();
    cout << "Overlap with Basis Set  : " << anm_list << " time: " << S.GetTime() << "(s)" << endl;
  }
  
  cout << endl << line << endl << endl;
  {
    StopWatch S;S.Start();
    const List<SMatrix<double> > anm_list  =db   ->GetRepulsion3C(*bs);
    S.Stop();
    cout << "Repulsion with Basis : " << anm_list << " time: " << S.GetTime() << "(s)" << endl;
  }
*/  
  
  
//---------------------------------------------------------------------------------------------------------------
//
//  2 and 3 center repulsion integrals
//
  if (repulsions)
  {
    cout << endl << line << endl << endl;
    ansm =db   ->GetRepulsion();
    StopWatch S;S.Start();
    numsm=dbnum->GetRepulsion();
    S.Stop();
    if (verbose) cout << "Analytic " << ansm << endl << "Numerical " << numsm << endl;
    cout << "Repulsion RMS error : " << RelError(ansm,numsm) << " time: " << S.GetTime() << "(s)" << endl;
  }
  if (repulsions)
  {
    cout << endl << line << endl << endl;
    anm =db   ->GetRepulsion();
    StopWatch S;S.Start();
    numm=dbnum->GetRepulsion(*bs);
    S.Stop();
    if (verbose) cout << "Analytic " << anm << endl << "Numerical " << numm << endl;
    cout << "Cross Repulsion RMS error : " << RelError(anm,numm) << " time: " << S.GetTime() << "(s)" << endl;
  }


  if (repulsions)
  {
    cout << endl << line << endl << endl;
    TBasisSetBrowser<double> bsbr(*bs);
    const ScalarFunction<double>& rsf(*bsbr);
    an =db   ->GetRepulsion(rsf);
    StopWatch S;S.Start();
    num=dbnum->GetRepulsion(rsf);
    S.Stop();
    if (verbose) cout << "Analytic " << an << endl << "Numerical " << num << endl;
    cout << "Repulsion with scalar Function RMS error : " << RelError(an,num) << " time: " << S.GetTime() << "(s)" << endl;
  }

		cout << endl << line << endl << endl;
  if (repulsions)
  {
    StopWatch S;S.Start();
    const MatrixList<double>& mlist  =db   ->GetRepulsion3C(*bs);
    S.Stop();
    cout << "Repulsion with Basis Set time: " << S.GetTime() << "(s)" << endl;
    S.Start();
    const MatrixList<double>& mlistnum  =dbnum   ->GetRepulsion3C(*bs);
    S.Stop();
    cout << "Numerial Repulsion with Basis Set time: " << S.GetTime() << "(s)" << endl;
    
    for (int i=0;i<bs->GetNumFunctions();i++)
    {
			if(verbose) cout << mlist[i] << mlistnum[i] << endl;
      cout << "Repulsion with Basis Set RMS error : " << RelError(mlist[i],mlistnum[i]) << endl;
    }
  }
  cout << endl << line << endl << endl;
//---------------------------------------------------------------------------------------------------------------
//
//  2 center nuclear and kinetic integrals
//
  cout << endl << line << endl << endl;
  {
    anm =db   ->GetNuclear(*cl);
    StopWatch S;S.Start();
    numm=dbnum->GetNuclear(*cl);
    S.Stop();
    if (verbose) cout << "Analytic " << anm << endl << "Numerical " << numm << endl;
    cout << "Nuclear Attraction RMS error : " << RelError(anm,numm) << " time: " << S.GetTime() << "(s)" << endl;
  }

  cout << endl << line << endl << endl;
  {
    anm =db   ->GetKinetic();
    StopWatch S;S.Start();
    numm=dbnum->GetKinetic();
    S.Stop();
    if (verbose) cout << "Analytic " << anm << endl << "Numerical " << numm << endl;
    cout << "Kinetic RMS error : " << RelError(anm,numm) << " time: " << S.GetTime() << "(s)" << endl;
  }


    
    

cleanup:
  delete bs;
  delete bsnum;
  delete cl;
    
  }
  return 0;
};


double RelError(const Vector<double>& a, const Vector<double>& b)
{
  return sqrt(SumSquares(a-b));
}

double RelError(const Matrix<double>& a, const Matrix<double>& b)
{
  
  return sqrt(SumSquares(a-b));
}

double RelError(const SMatrix<double>& a, const SMatrix<double>& b)
{
  
  return sqrt(SumSquares(a-b));
}
