#include <stdio.h>
#include <stdlib.h>
#include "atomic_functions.c"
#include "sle2.0.c"
#include "smplDynArray.c"
#include "B-splines.c"
#define KNODES 120



double f(double x)
{
	return -sin(x);
}

//double res[KNODES], Lsys[KNODES][KNODES], F[KNODES];
void solve_cup(double A, double B, int N)
{
	int i;
	double *res, *F,step = (B-A)/(N-1.0),arg, output;
	
	init1DArr(&res,N);
	init1DArr(&F,N);
	//init2DArr(&Lsys,N,N);
	
	for(i=0;i<N;i++)
		F[i] = f((double)(i)*step+A);
	
	printf("#Lsys init-ed\n");
	
	//solveLinearEquations(&Lsys, &F, &res, 0, N, 'D');
	solve3DiagInterp(f_dd_cup1(1.0),f_dd_cup1(0.0),&F,&res,N);
	
	printf("System solved\n");
	
	FILE *op;
	op = fopen("./output/plot", "w");
	for(arg=A;arg<B;arg+=0.01)
	{	
		//output = res[1]*f_cup((arg-A-step)/step)/1.8;
		output = 0.0;

		for(i=1;i<N-2;i++)
		{
			output += res[i]*f_cup(0.5*(arg-(double)(i)*step+A)/step)/3.714236763695158;
		}
		//output+= res[N-2]*f_cup((arg-B+step)/step)/1.8;
		fprintf(op,"%f %f\n",arg,output);
	}
	fclose(op);
	
	op = fopen("./output/dots", "w");
	for(arg=A;arg<=B;arg+=step)
	{
		fprintf(op,"%f %f\n",arg,sin(arg));
	}
	fclose(op);
	free1DArr(&res);
	free1DArr(&F);
}

void solve_B3(double A, double B, int N)
{
	int i;
	double *res, /* *Lsys,*/ *F,step = (B-A)/(N-1.0),arg, output;
	init1DArr(&res,N);
	init1DArr(&F,N);
	for(i=0;i<N;i++)
		F[i] = f((double)(i)*step+A);
		
	solve3DiagInterp(f_dd_B_3( 1.0),f_dd_B_3( 0.0),&F,&res,N);

	FILE *op;
	op = fopen("./output/plot", "w");
	for(arg=A;arg<B;arg+=0.01)
	{		
		output = -res[1]*f_B_3((arg-A)/step)*f_B_3(1.0)/6.0;
		//output = 0.0;
		for(i=1;i<=N-2;i++)
			output += res[i]*f_B_3((arg-(double)(i)*step+A)/step)/6.0;
		
		output+= -res[N-2]*f_B_3((arg-(double)(N-1)*step+A)/step)*f_B_3(1.0)/6.0;
		
		fprintf(op,"%f %f\n",arg,output);
	}
	fclose(op);
	
	op = fopen("./output/dots", "w");
	for(arg=A;arg<=B;arg+=step)
	{
		fprintf(op,"%f %f\n",arg,sin(arg));
	}
	fclose(op);
	free1DArr(&res);
	free1DArr(&F);
}

int main()
{
	init_up();
	init_cup();
	
	solve_B3(0.0,2*M_PI,KNODES);
	printf("%f\n",f_dd_cup1(0)/f_dd_cup2(0));
	/*
	double x=-2;
	for(x=-2.0;x<=2.0;x+=0.001)
		printf("%f %f %f %f %f\n",x,f_cup(x),f_d_cup(x),f_dd_cup1(x),f_dd_cup2(x));
	*/
	return 0;
}

