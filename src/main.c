#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
//#include "sle.c"
#include "sle2.0.c"
#include "smplDynArray.c"
#include <gsl/gsl_linalg.h>
#include "atomic_functions.c"
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



double d_f_e(double x)
{
	return cos(x);
}

double phi(double x, double y)
{
	return f_cup(x)*f_cup(y);
}

void insertA(gsl_matrix * sys, int i0, int j0, int n)
{
	int i;
	const double a=f_B_3(0.0),b=f_B_3(1.0);
	for(i=0;i<n;i++)
	{
		if(i==0)
		{
			gsl_matrix_set(sys, i0+0,j0+0,a);
			gsl_matrix_set(sys, i0+0,j0+1,b);
		}
		
		if(i!=0 && i!=n-1)
		{
			gsl_matrix_set(sys,i0+i,j0+i-1,b);
			gsl_matrix_set(sys,i0+i,j0+i,a);
			gsl_matrix_set(sys,i0+i,j0+i+1,b);
		}
		
		
		if(i==n-1)
		{
			gsl_matrix_set(sys,i0+n+1,j0+n-2,b);
			gsl_matrix_set(sys,i0+n+1,j0+n-1,a);
		}	
	}
}

void insertB(gsl_matrix * sys, int i0, int j0, int n)
{
	int i;
	const double a=f_B_3(1.0), b=f_B_3(sqrt(2.0));
	for(i=0;i<n;i++)
	{
		if(i==0)
		{
			gsl_matrix_set(sys, i0+0,j0+0,a);
			gsl_matrix_set(sys, i0+0,j0+1,b);
		}
		
		if(i!=0 && i!=n-1)
		{
			gsl_matrix_set(sys,i0+i,j0+i-1,b);
			gsl_matrix_set(sys,i0+i,j0+i,a);
			gsl_matrix_set(sys,i0+i,j0+i+1,b);
		}
		
		
		if(i==n-1)
		{
			gsl_matrix_set(sys,i0+n+1,j0+n-2,b);
			gsl_matrix_set(sys,i0+n+1,j0+n-1,a);
		}	
	}
}

void interp2D(double x1, double x2, double y1, double y2, int xs, int ys)
{
	int i,j,k;
	double outp, arg;
	const double 	xh=(x2-x1)/(xs-1.0), yh=(y2-y1)/(xs-1.0);
	//const double a=f_B_3(0.0),b=f_B_3(1.0),c=f_B_3(sqrt(2.0)) ;
	
	gsl_matrix * sys = gsl_matrix_alloc (xs*ys, xs*ys);
	gsl_vector * x = gsl_vector_alloc (xs*ys);
	gsl_vector * b = gsl_vector_alloc (xs*ys);
	
	for( i = 0; i < ys; i++ )
	{
		if(i==0)
		{
			insertA(sys,0,0,xs);
			insertB(sys,0,xs,xs);
		}		
		if(i!=0 && i!=ys-1)
		{
			insertB(sys,i*xs,(i-1)*xs,xs);
			insertA(sys,i*xs,i*xs,xs);
			insertB(sys,i*xs,(i+1)*xs,xs);
		}				
		if(i==ys-1)
		{
			insertB(sys,i*xs,(i-1)*xs,xs);
			insertA(sys,i*xs,i*xs,xs);
		}
	}
	
	
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
	gsl_vector_fprintf (stdout, x, "%g");
	
	FILE *op;
	op = fopen("./output/plot.b3", "w");
	/*for(i=0;i<n+2;i++)
		fprintf(op,"%f %f\n",(double)(i)*step,gsl_vector_get(x,i));
	*/
	
	for(arg = A; arg<=B; arg+=step)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_B_3((arg-step*(double)(i-1))/step);
			
		//fprintf(op,"%f %f\n",arg,fabs(outp-f_e(arg)));
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
		//fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	
	op = fopen("./output/plot.b3.d.err", "w");
	for(arg = A; arg<=B; arg+=0.01)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_d_B_3((arg-step*(double)(i-1))/step)/step;
			
		fprintf(op,"%f %f\n",arg,fabs(outp-d_f_e(arg)));
		//fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	op = fopen("./output/plot.b3.d", "w");
	for(arg = A; arg<=B; arg+=step)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_d_B_3((arg-step*(double)(i-1))/step)/step;
			
		//fprintf(op,"%f %f\n",arg,fabs(outp-d_f_e(arg)));
		fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	
	op = fopen("./output/plot.b3.dd.err", "w");
	for(arg = A; arg<=B; arg+=0.01)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_dd_B_3((arg-step*(double)(i-1))/step)/step/step;
			
		fprintf(op,"%f %f\n",arg,fabs(outp-f(arg)));
		//fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	op = fopen("./output/plot.b3.dd", "w");
	for(arg = A; arg<=B; arg+=0.01)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_dd_B_3((arg-step*(double)(i-1))/step)/step/step;
			
		//fprintf(op,"%f %f\n",arg,fabs(outp-d_f_e(arg)));
		fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	
	op = fopen("./output/solution", "w");
	for(arg=A;arg<=B;arg+=0.01)
	{
		fprintf(op,"%f %f\n",arg,f_e(arg));
	}
	fclose(op);
	op = fopen("./output/d_solution", "w");
	for(arg=A;arg<=B;arg+=0.01)
	{
		fprintf(op,"%f %f\n",arg,d_f_e(arg));
	}
	fclose(op);
	op = fopen("./output/dd_solution", "w");
	for(arg=A;arg<=B;arg+=0.01)
	{
		fprintf(op,"%f %f\n",arg,f(arg));
	}
	fclose(op);
	/*for(i=0;i<n+2;i++)
	{
		for(j=0;j<n+2;j++)
			printf("%5.2f ",sys[i][j]);
		printf("\t%5.2f\n", b[i]);
	}*/
	//solveLinearEquations(&sys, &b, &x, 0, n+2, 'D');
	//LUSolve(sys, b, &x, n+2);
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
			//b[0] = 0.0;
			gsl_matrix_set(sys, 0,0,f_cup(-1.0));
			gsl_matrix_set(sys, 0,1,f_cup(0.0));
			gsl_matrix_set(sys, 0,2,f_cup(1.0));
			
			//sys[0][1] = f_B_3(0.0);
			//sys[0][2] = f_B_3(1.0);
		}
		
		if(i!=0 && i!=n+1)
		{
			gsl_vector_set(b,i, f((double)(i-1)*step+A));
			gsl_matrix_set(sys,i,i-1,c1);
			gsl_matrix_set(sys,i,i,c0);
			gsl_matrix_set(sys,i,i+1,c1);
			
			//sys[i][i] 	= c0;
			//sys[i][i+1]	= c1;
		}
		
		
		if(i==n+1)
		{
			gsl_vector_set(b,i,0.0);
			gsl_matrix_set(sys,n+1,n-1,f_cup(-1.0));
			gsl_matrix_set(sys,n+1,n,f_cup(0.0));
			gsl_matrix_set(sys,n+1,n+1,f_cup(1.0));

			/*b[n+1] = 0.0;
			sys[n+1][n-1] = f_B_3(-1.0);
			sys[n+1][n] = f_B_3(0.0);
			sys[n+1][n+1] = f_B_3(1.0);*/
		}	
	}
	//gsl_matrix_fprintf (stdout, sys, "%g");
	gsl_permutation * p = gsl_permutation_alloc (n+2);
	gsl_linalg_LU_decomp (sys, p, &i);
	gsl_linalg_LU_solve (sys, p, b, x);
	gsl_vector_fprintf (stdout, x, "%g");
	
	FILE *op;
	op = fopen("./output/plot.cup", "w");
	/*for(i=0;i<n+2;i++)
		fprintf(op,"%f %f\n",(double)(i)*step,gsl_vector_get(x,i));
	*/
	
	for(arg = A; arg<=B; arg+=step)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_cup((arg-step*(double)(i-1))/step)/0.309520;
			
		//fprintf(op,"%f %f\n",arg,fabs(outp-f_e(arg)));
		fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	op = fopen("./output/plot.cup.d", "w");
	for(arg = A; arg<=B; arg+=step)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_d_cup((arg-step*(double)(i-1))/step)/sqrt(0.309520)/step;
			
		//fprintf(op,"%f %f\n",arg,fabs(outp-f_e(arg)));
		fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	op = fopen("./output/plot.cup.d.err", "w");
	for(arg = A; arg<=B; arg+=step)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_d_cup((arg-step*(double)(i-1))/step)/sqrt(0.309520)/step;
			
		fprintf(op,"%f %f\n",arg,fabs(outp-d_f_e(arg)));
		//fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	
	op = fopen("./output/plot.cup.dd", "w");
	for(arg = A; arg<=B; arg+=0.01)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_dd_cup1((arg-step*(double)(i-1))/step)/step/step;
			
		//fprintf(op,"%f %f\n",arg,fabs(outp-f_e(arg)));
		fprintf(op,"%f %f\n",arg,outp);
	}
	fclose(op);
	op = fopen("./output/plot.cup.dd.err", "w");
	for(arg = A; arg<=B; arg+=step)
	{
		outp = 0.0;
		for(i=0;i<n+2;i++)
			outp += gsl_vector_get(x,i)*f_dd_cup1((arg-step*(double)(i-1))/step)/0.309520/step/step*0.309519999;
			
		fprintf(op,"%f %f\n",arg,fabs(outp-f(arg)));
		//fprintf(op,"%f %f\n",arg,outp);
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
	/*for(i=0;i<n+2;i++)
	{
		for(j=0;j<n+2;j++)
			printf("%5.2f ",sys[i][j]);
		printf("\t%5.2f\n", b[i]);
	}*/
	//solveLinearEquations(&sys, &b, &x, 0, n+2, 'D');
	//LUSolve(sys, b, &x, n+2);
}

int main ()
{
	init_up();
	init_cup();
	
	solveB3(0.0,M_PI,1024);
	solveCup(0.0,M_PI,1024);
/*
	FILE *op;
	op = fopen("B3", "w");
	double x, y;
	for(x=-2.0; x<=2.0; x+=0.01)
		fprintf(op, "%f %f\n",x,f_d_B_3(x));
	fclose(op);	
	
	op = fopen("CUP", "w");
	//double x, y;
	for(x=-2.0; x<=2.0; x+=0.01)
		fprintf(op, "%f %f\n",x,f_d_cup(x));
	fclose(op);	*/
	return 0;
}
