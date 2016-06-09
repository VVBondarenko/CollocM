#include <stdio.h>
#include <math.h>
#define MaxProduct 20
/*	fourier-calculated atomic functions		*/
/*	version 0.1					*/
/*							*/
/*	structure: 					*/
/*	f_XXXX(double x, ... ) - atomic function "XXXX"	*/
/*	F_XXXX(double t, ... ) - image of "XXXX"*/
/*	o_XXXX(...) - any other function, used further	*/
/*							*/
/*							*/

double o_sinc(double x)
{
	if(x==0.)
		return 1.;
	else
		return sin(x)/x;
}


double F_up(double t)
{
	int k;
	double 	y=1., 
		coef = 2.;
	
	for(k=1;k<MaxProduct;k++)
	{
		y*=o_sinc(t/coef);
		coef*=2.;
	}
	return y;
}

double f_up(double x)
{
	if(fabs(x)>1.)
		return 0.;
	double 	pik,
		y = 0.5;
		
	for(pik = M_PI; pik<=M_PI*MaxProduct; pik+=M_PI)
	{
		y+=F_up(pik)*cos(pik*x);
	}
	return y;
}

double F_fup(double t, int n)
{
	return pow(o_sinc(t/2.),n)*F_up(t);
}

double f_fup(double x, int n)
{
	if((2.*fabs(x))>n+2)
	  	return 0.;
	  	
	double 	pik,
		y = 0.5,
		h = 2.*M_PI/((double)n +2.);
	for(pik = h; pik<=h*MaxProduct; pik+=h)
	{
		y+=F_fup(pik,n)*cos(pik*x);
	}
	return y*2./((double)n +2.);
}

double f_cup(double x)
{
	if(fabs(x)>2.)
		return 0.;
	double 	pik,
		y = 0.5;
		
	for(pik = M_PI; pik<=M_PI*MaxProduct; pik+=M_PI)
	{
		y+=pow(F_up(pik*0.5),2)*cos(pik*x*0.5);
	}
	return y*0.5;
}

double F_ha(double t, int a)
{
	int k;
	double 	y=1., 
		coef = a;
	
	for(k=1;k<MaxProduct;k++)
	{
		y*=o_sinc(t/coef);
		coef*=2.;
	}
	return y;
}

double f_ha(double x, int a)
{
	if(fabs(x)>1./(a-1))
		return 0.;
	
	double 	pik,
		y = 0.5,
		h = M_PI*(double)(a-1);
	for(pik = h; pik<=h*MaxProduct; pik+=h)
	{
		y+=F_ha(pik,a)*cos(pik*x);
	}
	return y*(double)(a-1);
}
/*
//used for debug
int main()
{
	FILE *graph;
	double i;
	graph = fopen("fup.dat", "w");
	for(i=-5.; i<=5.; i+=0.01)
		fprintf(graph,"%f %f\n", i, f_ha(i,5));
	fclose(graph);
	return 0;
}
*/
