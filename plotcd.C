// File: plotcd.C  plot a charge density.

#include "ChargeDensity/ChargeDensity.H"
#include "DFTDataBase/DFTDataBase.H"
#include "PlotterImplementation/PlotterImplementation.H"
#include "Mesh/LinearMesh.H"

#include <string>
#include <GetOpt.h>
#include <fstream.h>
#include <iostream.h>
#include "Misc/Unpickle.H"

string usage = " [-mMesh] ChargeDensity1 ChargeDensity2 ...\n";

int main(int argc, char **argv)
{
  if (argc==1) 
  {
    cerr << "usage: " << argv[0] << usage  << endl;
    return -1;
  }
  
//
//  Shared object database.
//
  string Dfile;
  DFTDataBase theDataBase(Dfile);
//
//  Intiate the charge density.
//
	string meshFile;
	SList<ChargeDensity> cds;
	Mesh*         mesh=0;
  StreamableObject::SetToBinary();
  
  GetOpt opt (argc, argv, "m:");
  int opchar;
  
  while ((opchar = opt ()) != EOF)
    switch (opchar)
    {
      case 'm': meshFile=opt.optarg; break;
      case '?': cerr << "usage: " << argv[0] << usage  << endl;return -1;
    }
//
//  Read in the mesh.
//
	if (!UnPickle(mesh,meshFile,"mesh") ) return -1;
	if (!mesh) mesh=new LinearMesh(-3,3,RVec(1,0,0),600);
//
//  Read in all the charge densities.
//
  for (int i=opt.optind;i<=argc-1;i++)
	{
		string cdFile(argv[i]);
		cout << "Reading charge density from file '" << cdFile << "'" << endl;
		ChargeDensity* cd;
		if (!UnPickle(cd  ,cdFile,"charge density") )
		{
			delete mesh;
			return -1;
		}
		cds.Add(cd);
	}
//  
//  Plot it.
//  

	{
		PlotterImplementation pl;
    pl.LogYAxis();
    uint n=1;
  
    pl.Symbols();
    pl.Plot(*cds[0],mesh,RVec(1,0,0));
//    pl.Lines();
//    pl.AddPlot(*cd);
    cin >> n;
  }
  
	delete mesh;
	
  return 0; 
}
