#include "pspline_lib.h"
#include "lecuyer.h"

struct L2_1D_DATA * dat_psgold;
struct SP_1D * sp_psgold;
double lambda_psgold;
int iknot_psgold;

double fminbr(double a,double b, double (*f)(double x),double tol);

double fknot_psgold(double t);
double flambda_psgold(double t);

double opt_knot_psgold();
double knot_loop_psgold(double shift, double scale);
double lambda_loop_psgold(double shift, double scale);

const int NUMBER_KNOT_ITER_PSGOLD = 10000;
const int NUMBER_KNOT_LOOPS_PSGOLD = 10;

const int NUMBER_LAMBDA_ITER_PSGOLD = 500;
const int NUMBER_LAMBDA_LOOPS_PSGOLD = 2;

const int MIN_SPACE_PSGOLD = 2;
const int NUMBER_LOOPS_PSGOLD = 20;

const double KNOT_TOL_PSGOLD = 0.001;
const double LAMBDA_TOL_PSGOLD = 0.001;

const double LAMBDA_MULT_PSGOLD = 0.9;
const double MIN_LOG_LAMBDA_PSGOLD = -23;
const double MAX_LOG_LAMBDA_PSGOLD = 5;

double slave_psgold(double shift, double scale, int seed, int stream);

void freepsgold(int* sampsize, double* xdata, double* ydata, int* ord,
      int* numknot, int* initseed, int* initstream, double* lambda,
      double* optknot, double* trace_hat, double* GCV, double* GSJS_value)
{
int i, order, n, nobs, inform, seed, stream;

double bestlambda;
double shift, scale;

struct L2_1D_DATA * origdat, * sortdat;

struct SP_1D * best_sp;

/*  set order and dimension for spline, seed, and stream */

order = *ord;
n = *numknot;
n += order;
seed = *initseed;
stream = *initstream;
  
/* load data into data struct */
origdat = data_1d_initialize(*sampsize, &inform);

for(i=0; i<*sampsize; i++)
{
    origdat->xdata[i] = xdata[i];
    origdat->ydata[i] = ydata[i];
}

nobs = origdat->nobs;

sortdat = data_1d_initialize(nobs, &inform);
dat_psgold = data_1d_initialize(nobs, &inform);
sort_data(origdat, sortdat);
scale_data(sortdat, dat_psgold, &shift, &scale);

sp_psgold = sp_1d_initialize( order, n,  &inform);
sp_1d_set_knots(sp_psgold, 0., 1.);
best_sp = sp_1d_initialize( order, n,  &inform);
sp_1d_set_knots(best_sp, 0., 1.);

bestlambda = slave_psgold(shift, scale, seed, stream);

/* Save information */
/* Lambda */
*lambda = bestlambda * pow(scale, 2*sp_psgold->order-2);

/* Optimal knots */ 
for (i = 0; i <*numknot; i++)
    optknot[i] = shift + scale * sp_psgold->knot[i+*ord];

/* trace of hat matrix */
*trace_hat = trace_hat_matrix(dat_psgold, sp_psgold, bestlambda);

/* GCV */
*GCV = gcv(dat_psgold, sp_psgold, bestlambda, MIN_SPACE_PSGOLD);

/* GSJS */
*GSJS_value = GSJS(dat_psgold);
}

double slave_psgold(double shift, double scale, int seed, int stream)
{
     int i;
     double currgcv, bestgcv, bestlambda = 0;
     
     struct SP_1D  * best_sp;       
     best_sp = spline_1d_copy(sp_psgold);

     RandomInitialise(stream, stream+1, seed);

     lambda_psgold = exp(RandomDouble(MIN_LOG_LAMBDA_PSGOLD,
                     MAX_LOG_LAMBDA_PSGOLD)); 
     bestgcv = INFINITY;
     sp_1d_set_knots(sp_psgold, 0., 1.);
     sp_1d_set_knots(best_sp, 0., 1.);

     for (i = 0; i < NUMBER_LOOPS_PSGOLD; i++)
     {
          currgcv = lambda_loop_psgold(shift, scale);
          if (currgcv < bestgcv)
          {
              bestlambda = lambda_psgold;
              bestgcv = currgcv;
          }
          currgcv = knot_loop_psgold(shift, scale);
          if (currgcv < bestgcv)
          {
              free_SP_1D(best_sp);
              best_sp = spline_1d_copy(sp_psgold);
              bestgcv = currgcv;
          }
     }

     return bestlambda;
}

double lambda_loop_psgold(double shift, double scale)
{
int i, j;

double glob_min_gcv, loc_min_gcv, currgcv;
double loc_best_lambda = 0, glob_best_lambda = 0;
double oldgcv, a, b;

glob_min_gcv = INFINITY;

for (j = 0; j < NUMBER_LAMBDA_LOOPS_PSGOLD; j++)
{
     loc_min_gcv = INFINITY;
     for (i = 0; i < NUMBER_LAMBDA_ITER_PSGOLD; i++)
     {
          lambda_psgold = exp(RandomDouble(MIN_LOG_LAMBDA_PSGOLD, MAX_LOG_LAMBDA_PSGOLD)); 
          currgcv = gcv_fit(dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD);
          if (currgcv < loc_min_gcv)
          {
              loc_best_lambda = lambda_psgold;
              loc_min_gcv = currgcv;
          }
      }

      lambda_psgold = loc_best_lambda;
      currgcv = loc_min_gcv;
      get_pen_1D_spline(dat_psgold, sp_psgold, lambda_psgold);

      if (currgcv < glob_min_gcv)
      {
          glob_best_lambda = lambda_psgold;
          glob_min_gcv = currgcv;
      }

      for(i=0; i<99; i++)
      {
          oldgcv = currgcv;
          a = log(lambda_psgold) + log(LAMBDA_MULT_PSGOLD);
          b = log(lambda_psgold) - log(LAMBDA_MULT_PSGOLD); 
          lambda_psgold = exp(fminbr(a, b, flambda_psgold, LAMBDA_TOL_PSGOLD));
          currgcv = gcv_fit( dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD );
          if( fabs( (currgcv-oldgcv)/(oldgcv + .0001) ) < .0001 )
	         break;
      }	

      currgcv = gcv_fit(dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD );

      if (currgcv < glob_min_gcv)
      {
          glob_best_lambda = lambda_psgold;
          glob_min_gcv = currgcv;
      }

}  /* end j loop */

lambda_psgold = glob_best_lambda;
currgcv = gcv_fit(dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD);

/* printf("The overall best GCV is %e\n", currgcv);
printf("The value of lambda is %e \n", lambda_psgold*pow(scale, 2*sp_psgold->order-2));
for(i=0; i< sp_psgold->n + sp_psgold->order; i++)
      printf("%f ", shift + scale * sp_psgold->knot[i]);
printf("\n\n"); */

return currgcv;
}


double knot_loop_psgold(double shift, double scale)
{
int i, j;

double glob_min_gcv, loc_min_gcv, currgcv;

struct SP_1D* loc_best_sp;
struct SP_1D* glob_best_sp;

loc_best_sp = spline_1d_copy(sp_psgold);
glob_best_sp = spline_1d_copy(sp_psgold);

glob_min_gcv = INFINITY;

for (j = 0; j < NUMBER_KNOT_LOOPS_PSGOLD; j++)
{
     loc_min_gcv = INFINITY;
     for (i = 0; i < NUMBER_KNOT_ITER_PSGOLD; i++)
     {
          /* generate random spline */
          get_random_knots(sp_psgold, dat_psgold);
     
          currgcv = gcv_fit(dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD);

          if (currgcv < loc_min_gcv)
          {
              free_SP_1D(loc_best_sp);
              loc_best_sp = spline_1d_copy(sp_psgold);
              loc_min_gcv = currgcv;
          }
      }
      
      free_SP_1D(sp_psgold);
      sp_psgold = spline_1d_copy(loc_best_sp);
      currgcv = loc_min_gcv;
      get_pen_1D_spline(dat_psgold, sp_psgold, lambda_psgold);

      if (currgcv < glob_min_gcv)
      {
          free_SP_1D(glob_best_sp);
          glob_best_sp = spline_1d_copy(sp_psgold);
          glob_min_gcv = currgcv;
      }

      currgcv = opt_knot_psgold();

      currgcv = gcv_fit(dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD );

/*    printf("The GCV is %e after opt \n", currgcv);
      printf("The value of lambda is %e \n", lambda_psgold*pow(scale, 2*sp_psgold->order-2));

      for(i=0; i< sp_psgold->n + sp_psgold->order; i++)
         printf("%f ", shift + scale * sp_psgold->knot[i]);
      printf("\n\n"); */

      if (currgcv < glob_min_gcv)
      {
          free_SP_1D(glob_best_sp);
          glob_best_sp = spline_1d_copy(sp_psgold);
          glob_min_gcv = currgcv;
      }

}  /* end j loop */

free_SP_1D(sp_psgold);
sp_psgold = spline_1d_copy(glob_best_sp);
currgcv = gcv_fit(dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD);

return currgcv;
}

double opt_knot_psgold(void){
	int   j;
	double oldgcv, currgcv, a, b; 
	
      oldgcv = gcv_fit( dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD);

      for(j=0; j<99; j++){
	for(iknot_psgold = sp_psgold->order; iknot_psgold < sp_psgold->n  ; iknot_psgold++){
		a = sp_psgold->knot[iknot_psgold-1], b = sp_psgold->knot[iknot_psgold+1];
		sp_psgold->knot[iknot_psgold] = fminbr(a, b, fknot_psgold, KNOT_TOL_PSGOLD);
	
	}

	for(iknot_psgold = sp_psgold->n -1; iknot_psgold >= sp_psgold->order  ; iknot_psgold--){
		a = sp_psgold->knot[iknot_psgold-1], b = sp_psgold->knot[iknot_psgold+1];
		sp_psgold->knot[iknot_psgold] = fminbr(a, b, fknot_psgold, KNOT_TOL_PSGOLD);
		
	}
	

	for(iknot_psgold = sp_psgold->n -1; iknot_psgold >= sp_psgold->order  ; iknot_psgold--){
		a = sp_psgold->knot[iknot_psgold-1], b = sp_psgold->knot[iknot_psgold+1];
		sp_psgold->knot[iknot_psgold] = fminbr(a, b, fknot_psgold, KNOT_TOL_PSGOLD);
		
	}

	for(iknot_psgold = sp_psgold->order; iknot_psgold < sp_psgold->n  ; iknot_psgold++){
		a = sp_psgold->knot[iknot_psgold-1], b = sp_psgold->knot[iknot_psgold+1];
		sp_psgold->knot[iknot_psgold] = fminbr(a, b, fknot_psgold, KNOT_TOL_PSGOLD);
	
	}
		currgcv = gcv_fit( dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD );
		if( fabs( (currgcv-oldgcv)/(oldgcv + .0001) ) < .0001 )
			break;
		oldgcv = currgcv;
}

	currgcv = gcv_fit( dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD );
	return currgcv;
}

double flambda_psgold(double t){
	double gcv_value;
      lambda_psgold = exp(t);
	gcv_value = gcv_fit(dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD);
	return gcv_value;
}

double fknot_psgold(double t){
	double currgcv;
	sp_psgold->knot[iknot_psgold] = t;
	currgcv = gcv_fit( dat_psgold, sp_psgold, lambda_psgold, MIN_SPACE_PSGOLD );
	return currgcv;
}

	
