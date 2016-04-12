#include <math.h>
#include <stdlib.h>

void solve3DiagInterp(double a, double b, double **F, double **X, int N)
{
	double 	//a = (*A)[0][1],
			//b = (*A)[0][0],
			c = a*b/(b-a)/(a+b);
	int i;
	
	printf("#way1 start\n");
	(*F)[0] /= b;
	for(i=1; i < N; i++)
		(*F)[i] = ((*F)[i]/a - (*F)[i-1])*c;
	printf("#way2 start\n");
	(*X)[N-1] = (*F)[N-1];
	for(i=N-2;i>=0;i--)
		(*X)[i] = (*F)[i] - c*(*X)[i+1];
	printf("#way2 end\n");
}
