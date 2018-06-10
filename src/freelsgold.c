#include "pspline_lib.h"
#include "lecuyer.h"

struct L2_1D_DATA * dat_lsgold;
struct SP_1D * sp_lsgold;
int iknot_lsgold;

void slave_lsgold(double shift, double scale, int seed, int stream);

double f_lsgold(double t);
double opt_knot_lsgold();
double knot_loop_lsgold(double shift, double scale);

const int NUMBER_ITER_LSGOLD = 10000;
const int NUMBER_LOOPS_LSGOLD = 200;
const int MIN_SPACE_LSGOLD = 2;

const double TOL_LSGOLD = 0.001;

double fminbr(double a,double b, double (*f)(double x),double tol);

void freelsgold(int* sampsize, double* xdata, double* ydata, int* ord,
      int* numknot, int* initseed, int* initstream, double* optknot, 
      double* trace_hat, double* GCV, double* GSJS_value)
{
int i, order, n, nobs, inform, seed, stream;

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
dat_lsgold = data_1d_initialize(nobs, &inform);
\
sort_data(origdat, sortdat);
scale_data(sortdat, dat_lsgold, &shift, &scale);

sp_lsgold = sp_1d_initialize( order, n,  &inform);
sp_1d_set_knots(sp_lsgold, 0., 1.);
best_sp = sp_1d_initialize( order, n,  &inform);
sp_1d_set_knots(best_sp, 0., 1.);

slave_lsgold(shift, scale, seed, stream);

/* Save information */
/* Optimal knots */ 
for (i = 0; i <*numknot; i++)
    optknot[i] = shift + scale * sp_lsgold->knot[i+*ord];

/* trace of hat matrix */
*trace_hat = trace_hat_matrix(dat_lsgold, sp_lsgold, 0.0);

/* GCV */
*GCV = gcv(dat_lsgold, sp_lsgold, 0.0, MIN_SPACE_LSGOLD);

/* GSJS */
*GSJS_value = GSJS(dat_lsgold);
}


void slave_lsgold(double shift, double scale, int seed, int stream)
{
     double currgcv, bestgcv;
 
     struct SP_1D  * best_sp;       
     best_sp = spline_1d_copy(sp_lsgold);

     RandomInitialise(stream, stream + 1, seed);
     
     sp_1d_set_knots(sp_lsgold, 0., 1.);
     sp_1d_set_knots(best_sp, 0., 1.);
     bestgcv = INFINITY;

     currgcv = knot_loop_lsgold(shift, scale);
     if (currgcv < bestgcv)
     {
         free_SP_1D(best_sp);
         best_sp = spline_1d_copy(sp_lsgold);
         bestgcv = currgcv;
     }
}

double knot_loop_lsgold(double shift, double scale)
{
int i, j;

double glob_min_gcv, loc_min_gcv, currgcv;

struct SP_1D* loc_best_sp;
struct SP_1D* glob_best_sp;

loc_best_sp = spline_1d_copy(sp_lsgold);
glob_best_sp = spline_1d_copy(sp_lsgold);

glob_min_gcv = INFINITY;

for (j = 0; j < NUMBER_LOOPS_LSGOLD; j++)
{
     loc_min_gcv = INFINITY;
     for (i = 0; i < NUMBER_ITER_LSGOLD; i++)
     {
          /* generate random spline */
          get_random_knots(sp_lsgold, dat_lsgold);

          /*solve L2 problem */
          get_L2_1D_spline( dat_lsgold, sp_lsgold );
     
          currgcv = gcv(dat_lsgold, sp_lsgold, 0.0, MIN_SPACE_LSGOLD);

          if (currgcv < loc_min_gcv)
          {
              free_SP_1D(loc_best_sp);
              loc_best_sp = spline_1d_copy(sp_lsgold);
              loc_min_gcv = currgcv;
          }
      }
      
      free_SP_1D(sp_lsgold);
      sp_lsgold = spline_1d_copy(loc_best_sp);
      currgcv = loc_min_gcv;
      get_L2_1D_spline(dat_lsgold, sp_lsgold);

      if (currgcv < glob_min_gcv)
      {
          free_SP_1D(glob_best_sp);
          glob_best_sp = spline_1d_copy(sp_lsgold);
          glob_min_gcv = currgcv;
      }

      currgcv = opt_knot_lsgold();

      get_L2_1D_spline(dat_lsgold, sp_lsgold);
      currgcv = gcv(dat_lsgold, sp_lsgold, 0.0, MIN_SPACE_LSGOLD );

      if (currgcv < glob_min_gcv)
      {
          free_SP_1D(glob_best_sp);
          glob_best_sp = spline_1d_copy(sp_lsgold);
          glob_min_gcv = currgcv;
      }

}  /* end j loop */

free_SP_1D(sp_lsgold);
sp_lsgold = spline_1d_copy(glob_best_sp);
get_L2_1D_spline(dat_lsgold, sp_lsgold);
currgcv = gcv(dat_lsgold, sp_lsgold, 0.0, MIN_SPACE_LSGOLD);

return currgcv;
}

double opt_knot_lsgold(void){
	int   j;
	double oldgcv, currgcv, a, b; 
	
      oldgcv = gcv( dat_lsgold, sp_lsgold, 0.0, MIN_SPACE_LSGOLD);

      for(j=0; j<99; j++){
	for(iknot_lsgold = sp_lsgold->order; iknot_lsgold < sp_lsgold->n  ; iknot_lsgold++){
		a = sp_lsgold->knot[iknot_lsgold-1], b = sp_lsgold->knot[iknot_lsgold+1];
		sp_lsgold->knot[iknot_lsgold] = fminbr(a, b, f_lsgold, TOL_LSGOLD);
	
	}

	for(iknot_lsgold = sp_lsgold->n -1; iknot_lsgold >= sp_lsgold->order  ; iknot_lsgold--){
		a = sp_lsgold->knot[iknot_lsgold-1], b = sp_lsgold->knot[iknot_lsgold+1];
		sp_lsgold->knot[iknot_lsgold] = fminbr(a, b, f_lsgold, TOL_LSGOLD);
		
	}
	

	for(iknot_lsgold = sp_lsgold->n -1; iknot_lsgold >= sp_lsgold->order  ; iknot_lsgold--){
		a = sp_lsgold->knot[iknot_lsgold-1], b = sp_lsgold->knot[iknot_lsgold+1];
		sp_lsgold->knot[iknot_lsgold] = fminbr(a, b, f_lsgold, TOL_LSGOLD);
		
	}

	for(iknot_lsgold = sp_lsgold->order; iknot_lsgold < sp_lsgold->n  ; iknot_lsgold++){
		a = sp_lsgold->knot[iknot_lsgold-1], b = sp_lsgold->knot[iknot_lsgold+1];
		sp_lsgold->knot[iknot_lsgold] = fminbr(a, b, f_lsgold, TOL_LSGOLD);
	
	}
		currgcv = gcv( dat_lsgold, sp_lsgold, 0.0, MIN_SPACE_LSGOLD );
		if( fabs( (currgcv-oldgcv)/(oldgcv + .0001) ) < .0001 )
			break;
		oldgcv = currgcv;
}

	currgcv = gcv( dat_lsgold, sp_lsgold, 0.0, MIN_SPACE_LSGOLD );
	return currgcv;
}

double f_lsgold(double t){
	double currgcv;
	sp_lsgold->knot[iknot_lsgold] = t;
	get_L2_1D_spline(dat_lsgold, sp_lsgold);
	currgcv = gcv( dat_lsgold, sp_lsgold, 0.0, MIN_SPACE_LSGOLD );
	return currgcv;
}
