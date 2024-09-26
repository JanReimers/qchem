#include <iostream.h>
#include <math.h>

int main()
{
	double del=0.0001,I=0;
	for (int n=0;n<=10/del;n++)
	{
		double r=n*del;
		if (n%100==1) cout << r << " " << I << " " << I/r << " " << I*r << endl;
		I=I*exp( (n*n-(n+1)*(n+1)) *del*del )+del;
	}
	return -1;
}
