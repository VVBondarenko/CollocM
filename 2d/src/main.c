#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>


#include "B-splines.c"
#include "af_fourier.c"
int         N = 5;
double      X0 = -1.,
			X1 = 1.,
			Y0 = -1.,
			Y1 = 1.,
			stepx, stepy, diff_step, glob_delta;


double right_part_f(double x, double y)
{
	//return 12.*(y*y*(x*x*x*x-1.) + x*x*(y*y*y*y-1.));
	//return -2.*sin(x)*sin(y);
	return 2.*exp((-1.+x*x)*(-1.+y*y))*(-2.+3.*y*y+2.*x*x*x*x*y*y + 
     x*x*(3. - 8.*y*y + 2.*y*y*y*y));
}

double u_exact(double x, double y)
{
	return exp((x*x-1.)*(y*y-1.))-1.;
	
	//return sin(x)*sin(y);
}

double omega(double x, double y)
{
	//double h = (x-X0)*(x-X1)*(y-Y0)*(y-Y1);	
	//if(h>=1.) return 1.;
	//else
	//{
		//return h;
	//}
	return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
	//return (1-x)*(1+x)+(1-y)*(1+y)-sqrt(pow(1-x,2)*pow(1+x,2)+pow(1-y,2)*pow(1+y,2));
	//return tanh((x-X0)*(x-X1)*(y-Y0)*(y-Y1)*10.);
	//return tanh(100.*(x-X0))*tanh(100.*(X1-x))*tanh(100.*(y-Y0))*tanh(100.*(Y1-y));
	//return (x*x-1.)*(y*y-1.);
	//return 1.0;
}

double phi(double x, double y)
{
	//return f_B_3(x)*f_B_3(y);
	//return f_fup(x,2)*f_fup(y,2);
	//return f_cup(x)*f_cup(y);
	return f_fup3_poly(x)*f_fup3_poly(y);
}

double phi_fup(double x, double y)
{
	//return f_B_3(x)*f_B_3(y);
	return f_fup(x,3)*f_fup(y,3);
	//return f_cup(x)*f_cup(y);
}

double basis(double x, double y, int n)
{
/*
 * 
 * 	- \phi*\omega везде					fup	0.025742 2/3 coef | 1/3 0.040864
 *	- \phi*\omega только на границе		fup	0.025629 2/3 coef | 1/3 0.137500

 *	- на границе						B_3	0.294823 2/3 coef | 1/3 0.022808
 *	- везде								B_3	0.293138 2/3 coef | 1/3 0.007089
 * 
 */
	//if( n/N <= 2 || n%N <= 2 || n%N >= (N-3) || n/N >= (N-3))
		return phi( 0.33333333333*(x-X0-stepx*(double)(1+n/N))/stepx,
			0.33333333333*(y-Y0-stepy*(double)(1+n%N))/stepy )* omega(x,y);
	//else
		//return phi( 0.33333333333*(x-X0-stepx*(double)(1+n/N))/stepx,
			//0.33333333333*(y-Y0-stepy*(double)(1+n%N))/stepy );
			
}

double laplacian_basis(double x, double y, int n)
{
	return glob_delta*glob_delta*(basis(x+diff_step,y,n)+basis(x-diff_step,y,n)+
						basis(x,y+diff_step,n)+basis(x,y-diff_step,n)
					  -4.*basis(x,y,n));
}

void form_matrix
     (gsl_matrix * system,
	gsl_vector * RightPart,
	double x1, double x2,
	double y1, double y2)
{
	int i, j;
	for(i = 0; i < N*N; i++)
	{
		gsl_vector_set(RightPart, i, right_part_f(x1+stepx*(double)(1+i/N),
									y1+stepy*(double)(i%N+1)));
		for(j = 0; j < N*N; j++)
		{
			gsl_matrix_set(system, i,j, laplacian_basis(
											x1+stepx*(double)(1+i/N),
											y1+stepy*(double)(i%N+1), j));
		}
	}
}

void solve_matrix_eq
     (gsl_vector * solution,
	gsl_matrix * system,
	gsl_vector * RightPart)
//Solve SLE Ax=b, where A = system, b = RightPart, x = solution
{
	int i;

	gsl_permutation * p = gsl_permutation_alloc (N*N);
	gsl_linalg_LU_decomp (system, p, &i);
	gsl_linalg_LU_solve (system, p, RightPart, solution);
}

double reconstruct_at(gsl_vector *solution, double x, double y)
// Reconstucts value of solution at point (x,y)
{
	int i; double result = 0.;
	for(i=0; i<N*N; i++)
	{
		result += gsl_vector_get(solution, i)*basis(x,y,i);
		//result = fmax(basis(x,y,i),result);
	}
	return result;
}

void plot_region
     (gsl_vector *solution,
	double x1, double x2,
	double y1, double y2)
// Plot solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y

{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j;
	FILE * op;
	op = fopen("../plot_data/plot_region", "w");
	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			fprintf(op, "%f %f %f\n", i,j, reconstruct_at(solution,i,j));
	fclose(op);
	i = system("./Plot");
}

void plot_region_error
     (gsl_vector *solution, 
      double x1, double x2, 
      double y1, double y2)
// Plot abs error of solution in rectangle region 
// from x1 till x2 by x, and from y1 till y2 by y
{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j;
	FILE * op;
	op = fopen("../plot_data/plot_region_error", "w");
	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			fprintf(op, "%f %f %f\n", i,j, fabs(reconstruct_at(solution,i,j)-u_exact(i,j)));
	fclose(op);
	i = system("./Plot_err");
}

void errors_to_stdio
     (gsl_vector *solution, 
      double x1, double x2, 
      double y1, double y2)
{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j,maxerr=0.,err;

	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			if((err = fabs(reconstruct_at(solution,i,j)-u_exact(i,j)))>maxerr) maxerr=err;
	printf("%d %f\n",N*N,maxerr);
}

int main(int argc, char **argv)
{
	if(argc > 0)
		N = atoi(argv[1]);
	stepx = (X1-X0)/(double)(2*N);
	stepy = (Y1-Y0)/(double)(2*N);
	N = 2*N-1;
	diff_step = pow(2.,-3);
	glob_delta = 1./diff_step;


	gsl_matrix 	*sys 		= gsl_matrix_alloc (N*N,N*N);;
	gsl_vector  *rightpart	= gsl_vector_alloc(N*N),
			*solution	= gsl_vector_alloc(N*N);

	form_matrix			(sys, rightpart, X0,X1, Y0,Y1);
	//printf("formed\n");
	//FILE *op;
	//op = fopen("./matrix", "w");
	//gsl_matrix_fprintf	(op, sys, "%f");
	//fclose(op);
	
	solve_matrix_eq		(solution, sys, rightpart);
	
	plot_region			(solution, X0,X1, Y0,Y1);
	plot_region_error	(solution, X0,X1, Y0,Y1);
	
	errors_to_stdio		(solution, X0,X1, Y0,Y1);
	
	return 0;
}
