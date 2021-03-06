#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
//#include "sle.c"
#include "sle2.0.c"
#include "smplDynArray.c"
#include <gsl/gsl_linalg.h>
#include "atomic_functions.c"
#include <string.h>


/*
double f(double x)
{
	return -sin(x);
}

double f_e(double x)
{
	return sin(x);
}
*/
double f(double x)
{
	return -exp(x*(x-M_PI))*((2.0*x-M_PI)*(2.0*x-M_PI)+2.0);
}

double f_e(double x)
{
	return 1.0 - exp(x*(x-M_PI));
}



void solveB3(double A, double B, int n)
{
	int i,j;
	double outp, arg;
	const double step1=(B-A)/(n+1.0), step=(B-A)/(n-1.0), c0=f_dd_B_3(0.0)/step/step,c1=f_dd_B_3(1.0)/step/step ;
	
	gsl_matrix * sys = gsl_matrix_alloc (n+2, n+2);
	gsl_vector * x = gsl_vector_alloc (n+2);
	gsl_vector * b = gsl_vector_alloc (n+2);

	for(i=0;i<n+2;i++)
	{
		if(i==0)
		{
			gsl_vector_set(b,0,0.0);
			gsl_matrix_set(sys, 0,0,f_B_3(-1.0));
			gsl_matrix_set(sys, 0,1,f_B_3(0.0));
			gsl_matrix_set(sys, 0,2,f_B_3(1.0));

		}
		
		if(i!=0 && i!=n+1)
		{
			gsl_vector_set(b,i, f((double)(i-1)*step+A));
			gsl_matrix_set(sys,i,i-1,c1);
			gsl_matrix_set(sys,i,i,c0);
			gsl_matrix_set(sys,i,i+1,c1);
		}
		
		
		if(i==n+1)
		{
			gsl_vector_set(b,i,0.0);
			gsl_matrix_set(sys,n+1,n-1,f_B_3(-1.0));
			gsl_matrix_set(sys,n+1,n,f_B_3(0.0));
			gsl_matrix_set(sys,n+1,n+1,f_B_3(1.0));
		}	
	}
	//gsl_matrix_fprintf (stdout, sys, "%g");
	gsl_permutation * p = gsl_permutation_alloc (n+2);
	gsl_linalg_LU_decomp (sys, p, &i);
	gsl_linalg_LU_solve (sys, p, b, x);
	//gsl_vector_fprintf (stdout, x, "%g");
	
	FILE *op;
	
	op = fopen("./output/plot.b3", "w");
	for(arg = A; arg<=B; arg+=step)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_B_3((arg-step*(double)(i-1))/step);
		fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	op = fopen("./output/plot.b3.err", "w");
	for(arg = A; arg<=B; arg+=0.01)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_B_3((arg-step*(double)(i-1))/step);
			
		fprintf(op,"%f %f\n",arg,fabs(outp-f_e(arg)));
	}
	fclose(op);

	
	op = fopen("./output/solution", "w");
	for(arg=A;arg<=B;arg+=0.01)
	{
		fprintf(op,"%f %f\n",arg,f_e(arg));
	}
	fclose(op);
	
}

void solveCup(double A, double B, int n)
{
	int i,j;
	double outp, arg;
	const double step1=(B-A)/(n+1.0), step=(B-A)/(n-1.0), c0=f_dd_cup1(0.0)/step/step,c1=f_dd_cup1(1.0)/step/step ;
	
	gsl_matrix * sys = gsl_matrix_alloc (n+2, n+2);
	gsl_vector * x = gsl_vector_alloc (n+2);
	gsl_vector * b = gsl_vector_alloc (n+2);

	for(i=0;i<n+2;i++)
	{
		if(i==0)
		{
			gsl_vector_set(b,0,0.0);
			gsl_matrix_set(sys, 0,0,f_cup(-1.0));
			gsl_matrix_set(sys, 0,1,f_cup(0.0));
			gsl_matrix_set(sys, 0,2,f_cup(1.0));
		}
		
		if(i!=0 && i!=n+1)
		{
			gsl_vector_set(b,i, f((double)(i-1)*step+A));
			gsl_matrix_set(sys,i,i-1,c1);
			gsl_matrix_set(sys,i,i,c0);
			gsl_matrix_set(sys,i,i+1,c1);
		}
		
		
		if(i==n+1)
		{
			gsl_vector_set(b,i,0.0);
			gsl_matrix_set(sys,n+1,n-1,f_cup(-1.0));
			gsl_matrix_set(sys,n+1,n,f_cup(0.0));
			gsl_matrix_set(sys,n+1,n+1,f_cup(1.0));
		}	
	}
	//gsl_matrix_fprintf (stdout, sys, "%g");
	gsl_permutation * p = gsl_permutation_alloc (n+2);
	gsl_linalg_LU_decomp (sys, p, &i);
	gsl_linalg_LU_solve (sys, p, b, x);
	//gsl_vector_fprintf (stdout, x, "%g");
	
	FILE *op;
	
	op = fopen("./output/plot.cup", "w");
	for(arg = A; arg<=B; arg+=step)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_cup((arg-step*(double)(i-1))/step)/0.309520;
			
		//fprintf(op,"%f %f\n",arg,fabs(outp-f_e(arg)));
		fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	

	op = fopen("./output/plot.cup.err", "w");
	for(arg = A; arg<=B; arg+=0.01)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_cup((arg-step*(double)(i-1))/step)/0.309520;
			
		fprintf(op,"%f %f\n",arg,fabs(outp-f_e(arg)));
		//fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	
	op = fopen("./output/dots", "w");
	for(arg=A;arg<=B;arg+=0.01)
	{
		fprintf(op,"%f %f\n",arg,f_e(arg));
	}
	fclose(op);

}


int main (int arc, char** argv)
{
	init_up();
	init_cup();
	
	
	int resol = atoi(argv[1]);
	solveB3(0.0,M_PI,resol);
	solveCup(0.0,M_PI,resol);
	
	//plotting error graphics
	char args[80];
	strcpy(args,"gnuplot -e \"outputname='./img/err");
	strcat(args,argv[1]);
	strcat(args,".eps'\" errplot2f.conf");
	system(args);
	
	//plotting graphics
	strcpy(args,"gnuplot -e \"outputname='./img/plot");
	strcat(args,argv[1]);
	strcat(args,".eps'\" plot2f.conf");
	system(args);
	printf("%f %f %f\n",f_B_3(0.),f_B_3(1.),f_B_3(2.));
	return 0;
}
