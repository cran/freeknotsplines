#include "pspline_lib.h"
#include "lecuyer.h"

struct L2_1D_DATA *  dat_psgen;
struct SP_1D * sp_psgen;
double lambda_psgen;

const int KNOT_POP_SIZE_PSGEN = 1000;
const int KNOT_NUMBER_FIT_PSGEN = 100;
const int KNOT_NUMBER_MUTATE_PSGEN = 1000;
const int KNOT_NUMBER_ITER_PSGEN = 200;

const int LAMBDA_POP_SIZE_PSGEN = 100;
const int LAMBDA_NUMBER_FIT_PSGEN = 10;
const int LAMBDA_NUMBER_MUTATE_PSGEN = 100;
const int LAMBDA_NUMBER_ITER_PSGEN = 100;

const double LAMBDA_CHILD_MULT_PSGEN = 0.9;
const double MIN_LOG_LAMBDA_PSGEN = -23;
const double MAX_LOG_LAMBDA_PSGEN = 5;

const int MIN_SPACE_PSGEN = 2;
const int NUMBER_LOOPS_PSGEN = 20;

const double END_EPSILON_PSGEN = 0.001;
const int END_CHECK_PSGEN = 10;

void lambda_reproduce_psgen(double * lambda1, double * lambda2, double * child);
void lambda_mutate_psgen(double * child);

void knot_mutate_psgen(struct SP_1D * child);
void knot_reproduce_psgen(struct SP_1D * par1, struct SP_1D * par2, struct SP_1D * child);

double lambda_loop_psgen(double shift, double scale);
double knot_loop_psgen(double shift, double scale);

int qsort_knot_partition_psgen(double * gcv_array, struct SP_1D ** par_array, 
         int left, int right, int pivotIndex);
void qsort_knot_psgen(double * gcv_array, struct SP_1D ** par_array, 
         int left, int right, int numSort);

int qsort_lambda_partition_psgen(double * gcv_array, double * par_array, 
          int left, int right, int pivotIndex);
void qsort_lambda_psgen(double * gcv_array, double * par_array, int left, int right, int numSort);

double slave_psgen(double shift, double scale, int seed, int stream);

void freepsgen(int* sampsize, double* xdata, double* ydata, int* ord,
      int* numknot, int* initseed, int* initstream, double* lambda,
      double* optknot, double* trace_hat, double* GCV, double* GSJS_value)
{
int i, order, n, nobs, inform, seed, stream;

double shift, scale;
double bestlambda;

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
dat_psgen = data_1d_initialize(nobs, &inform);
sort_data(origdat, sortdat);
scale_data(sortdat, dat_psgen, &shift, &scale);

sp_psgen = sp_1d_initialize( order, n,  &inform);
sp_1d_set_knots(sp_psgen, 0., 1.);
best_sp = sp_1d_initialize( order, n,  &inform);
sp_1d_set_knots(best_sp, 0., 1.);

bestlambda = slave_psgen(shift, scale, seed, stream);

/* Save information */
/* Lambda */
*lambda = bestlambda * pow(scale, 2*sp_psgen->order-2);

/* Optimal knots */ 
for (i = 0; i <*numknot; i++)
    optknot[i] = shift + scale * sp_psgen->knot[i+*ord];

/* trace of hat matrix */
*trace_hat = trace_hat_matrix(dat_psgen, sp_psgen, bestlambda);

/* GCV */
*GCV = gcv(dat_psgen, sp_psgen, bestlambda, MIN_SPACE_PSGEN);

/* GSJS */
*GSJS_value = GSJS(dat_psgen);
}

double slave_psgen(double shift, double scale, int seed, int stream)
{
     int i;
     double currgcv, bestgcv, bestlambda = 0;
  
     struct SP_1D  * best_sp;       
     best_sp = spline_1d_copy(sp_psgen);
  
     RandomInitialise(stream, stream + 1, seed);
     sp_1d_set_knots(sp_psgen, 0., 1.);
     sp_1d_set_knots(best_sp, 0., 1.);
     bestgcv = INFINITY;

     for (i = 0; i < NUMBER_LOOPS_PSGEN; i++)
     {
          currgcv = lambda_loop_psgen(shift, scale);
          if (currgcv < bestgcv)
          {
              bestlambda = lambda_psgen;
              bestgcv = currgcv;
          }
          currgcv = knot_loop_psgen(shift, scale);
          if (currgcv < bestgcv)
          {
              free_SP_1D(best_sp);
              best_sp = spline_1d_copy(sp_psgen);
              bestgcv = currgcv;
          }
     }

     return bestlambda;
}

double lambda_loop_psgen(double shift, double scale)
{
int i, j;
int rand1, rand2;

double bestgcv;
double best_lambda = 0;

double par_array[LAMBDA_POP_SIZE_PSGEN];
double child_array[LAMBDA_POP_SIZE_PSGEN];
double gcv_array[LAMBDA_POP_SIZE_PSGEN];

double old_crit, new_crit, conv_crit;
int end_flag = 0;

bestgcv = INFINITY;

old_crit = INFINITY;

/* generate random lambdas for parents */
for (i = 0; i < LAMBDA_POP_SIZE_PSGEN; i++)
{
    par_array[i] = exp(RandomDouble(MIN_LOG_LAMBDA_PSGEN, MAX_LOG_LAMBDA_PSGEN)); 
}

/* Begin major loop  */
for(j=0; j< LAMBDA_NUMBER_ITER_PSGEN; j++){

/* solve L2 problem for parents and compute errors for parents */
for (i = 0; i < LAMBDA_POP_SIZE_PSGEN; i++)
{
    gcv_array[i] = gcv_fit( dat_psgen, sp_psgen, par_array[i], MIN_SPACE_PSGEN);
}

/* sort parent splines by gcv */
/* note: only the fittest splines need to be sorted */
qsort_lambda_psgen(gcv_array, par_array, 0, LAMBDA_POP_SIZE_PSGEN-1, LAMBDA_NUMBER_FIT_PSGEN);

/* update best spline */
if(gcv_array[0] < bestgcv){
	bestgcv = gcv_array[0];
	best_lambda = par_array[0];
}

/* update termination criterion */
if (j % END_CHECK_PSGEN == 0)
{
     lambda_psgen = par_array[0];
     get_pen_1D_spline(dat_psgen, sp_psgen, lambda_psgen);
     new_crit = gcv_array[0] - get_L2_error(sp_psgen, dat_psgen)/(dat_psgen->nobs);
     conv_crit  = fabs(new_crit - old_crit)/old_crit;
     if (conv_crit < END_EPSILON_PSGEN)
          end_flag++;  /* if equal to 0, raise to 1; if equal to 1, raise to 2 */
     else
          if (end_flag == 1) /* reset end_flag to 0 */
             end_flag = 0;
     old_crit = new_crit;
}

/* terminate if necessary */
if (end_flag == 2)
{
    break;
}

/* create child splines */
for (i = 0; i < LAMBDA_POP_SIZE_PSGEN; i++)
{
     rand1 = RandomInt(0, LAMBDA_NUMBER_FIT_PSGEN - 1);
     rand2 = RandomInt(0, LAMBDA_NUMBER_FIT_PSGEN - 1);
     lambda_reproduce_psgen(&par_array[rand1], &par_array[rand2], &child_array[i]);
}

/* mutate child splines */
for (i = 0; i < LAMBDA_NUMBER_MUTATE_PSGEN; i++)
{
     rand1 = RandomInt(0, LAMBDA_POP_SIZE_PSGEN - 1);
     lambda_mutate_psgen(&child_array[rand1]);
}

/* copy child splines to parent splines */
for (i = 0; i < LAMBDA_POP_SIZE_PSGEN; i++)
    par_array[i] = child_array[i];

} /* j loop completed*/

lambda_psgen = best_lambda;
get_pen_1D_spline(dat_psgen, sp_psgen, lambda_psgen);

return bestgcv;
}

double knot_loop_psgen(double shift, double scale)
{
int i,j;
int inform;
int rand1, rand2;

struct SP_1D  * best_sp;
struct SP_1D  * par_array[KNOT_POP_SIZE_PSGEN];
struct SP_1D  * child_array[KNOT_POP_SIZE_PSGEN];

double  gcv_array[KNOT_POP_SIZE_PSGEN];
double  bestgcv; 

double old_crit, new_crit, conv_crit;
int end_flag = 0;

best_sp = spline_1d_copy(sp_psgen);
bestgcv = INFINITY;

old_crit = INFINITY;

/* initialize parent splines */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
    par_array[i] = sp_1d_initialize( sp_psgen->order, sp_psgen->n,  &inform);

for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
    sp_1d_set_knots(par_array[i], 0., 1.);

/* generate random splines for parents */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
    get_random_knots(par_array[i], dat_psgen);

/* Begin major loop  */
for(j=0; j< KNOT_NUMBER_ITER_PSGEN; j++){

/* solve L2 problem for parents and compute errors for parents */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
{
    gcv_array[i] = gcv_fit( dat_psgen, par_array[i], lambda_psgen, MIN_SPACE_PSGEN);
}

/* sort parent splines by psse */
/* note: only the fittest splines need to be sorted */
qsort_knot_psgen(gcv_array, par_array, 0, KNOT_POP_SIZE_PSGEN-1, KNOT_NUMBER_FIT_PSGEN);

/* update best spline */
if(gcv_array[0] < bestgcv){
	free_SP_1D(best_sp);
	bestgcv = gcv_array[0];
	best_sp = spline_1d_copy(par_array[0]);
}

/* update termination criterion */
if (j % END_CHECK_PSGEN == 0)
{
     new_crit = gcv_array[0] - get_L2_error(best_sp, dat_psgen)/(dat_psgen->nobs);
     conv_crit  = fabs(new_crit - old_crit)/old_crit;
     if (conv_crit < END_EPSILON_PSGEN)
          end_flag++;  /* if equal to 0, raise to 1; if equal to 1, raise to 2 */
     else
          if (end_flag == 1) /* reset end_flag to 0 */
             end_flag = 0;
     old_crit = new_crit;
}

/* terminate if necessary */
if (end_flag == 2)
{
    break;
}

/* initialize child splines */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
     child_array[i] = sp_1d_initialize( sp_psgen->order, sp_psgen->n,  &inform);

for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
     sp_1d_set_knots(child_array[i], 0., 1.);

/* create child splines */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
{
     rand1 = RandomInt(0, KNOT_NUMBER_FIT_PSGEN - 1);
     rand2 = RandomInt(0, KNOT_NUMBER_FIT_PSGEN - 1);
     knot_reproduce_psgen(par_array[rand1], par_array[rand2], child_array[i]);
}

/* mutate child splines */
for (i = 0; i < KNOT_NUMBER_MUTATE_PSGEN; i++)
{
     rand1 = RandomInt(0, KNOT_POP_SIZE_PSGEN - 1);
     knot_mutate_psgen(child_array[rand1]);
}

/* free parent splines */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
    free_SP_1D(par_array[i]);

/* copy child splines to parent splines */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
    par_array[i] = spline_1d_copy(child_array[i]);

/* free child splines */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
    free_SP_1D(child_array[i]);
}

/* j loop completed*/
/* free parent splines */
for (i = 0; i < KNOT_POP_SIZE_PSGEN; i++)
    free_SP_1D(par_array[i]);

free_SP_1D(sp_psgen);
sp_psgen = spline_1d_copy(best_sp);

return bestgcv;
}

void lambda_reproduce_psgen(double *lambda1, double * lambda2, double * child)
{
      double max, min;
      max = (*lambda1 > *lambda2)? *lambda1 : *lambda2;
      min = *lambda1 + *lambda2 - max;
      *child = exp(RandomDouble(log(min), log(max)));
}

void lambda_mutate_psgen(double * child)
{
     double unif = RandomDouble(log(LAMBDA_CHILD_MULT_PSGEN), -log(LAMBDA_CHILD_MULT_PSGEN));
     *child *= exp(unif); 
}

void knot_reproduce_psgen(struct SP_1D * par1, struct SP_1D * par2, struct SP_1D * child)
{
     int i, knot_break, knot_index1, knot_index2;
     child->n = par1->n;
     child->order = par1->order;
     for (i = 0; i < child->order; i++)
          child->knot[i] = 0.0;

      /* first knot in par2 */
     knot_break = RandomInt(child->order, child->n);
     knot_index1 = child->order;
     knot_index2 = knot_break;
     
     for (i = child->order; i < child->n; i++)
     {
          if ((knot_index1 == knot_break) ||
                 (par2->knot[knot_index2] <= par1->knot[knot_index1]))
          {
              child->knot[i] = par2->knot[knot_index2];
              knot_index2++;
          }
          else if ((knot_index2 == child->n) ||
                   (par1->knot[knot_index1] < par2->knot[knot_index2]))
          {
              child->knot[i] = par1->knot[knot_index1];
              knot_index1++;
          }
     }
     for (i = child->n; i < child->n+child->order; i++)
          child->knot[i] = 1.0;
}

void knot_mutate_psgen(struct SP_1D * child)
{
     int knot_number;
     knot_number = RandomInt(child->order, child->n - 1);
     child->knot[knot_number] = RandomDouble(child->knot[knot_number-1], 
           child->knot[knot_number+1]);
}

int qsort_lambda_partition_psgen(double * gcv_array, double * par_array, int left, int right, int pivotIndex)
{
     int storeIndex, i;
     double pivotValue;
     
     pivotValue = gcv_array[pivotIndex];

     /* Move pivot to end */
     swap_scalar(&(gcv_array[pivotIndex]), &(gcv_array[right]));
     swap_scalar(&(par_array[pivotIndex]), &(par_array[right]));

     storeIndex = left - 1;

     for(i = left; i <= right-1; i++)
         if (gcv_array[i] <= pivotValue)
         {
             storeIndex++;
             swap_scalar(&(gcv_array[storeIndex]), &(gcv_array[i]));
             swap_scalar(&(par_array[storeIndex]), &(par_array[i]));
         }

     /* Move pivot to its final place */
     swap_scalar(&(gcv_array[right]), &(gcv_array[storeIndex + 1]));
     swap_scalar(&(par_array[right]), &(par_array[storeIndex + 1]));

     return (storeIndex + 1);
}

void qsort_lambda_psgen(double * gcv_array, double * par_array, int left, int right, int numSort)
{
     /* Only sorts first numSort items */

     int pivotIndex, pivotNewIndex;

     if (right > left)
     {
           pivotIndex = left;
           pivotNewIndex = qsort_lambda_partition_psgen(gcv_array, par_array, 
                         left, right, pivotIndex);
           qsort_lambda_psgen(gcv_array, par_array, left, pivotNewIndex - 1, numSort);
           if (pivotNewIndex < numSort)
                 qsort_lambda_psgen(gcv_array, par_array, pivotNewIndex + 1, right, numSort);
     }
}

int qsort_knot_partition_psgen(double * gcv_array, struct SP_1D ** par_array, int left, int right, int pivotIndex)
{
     int storeIndex, i;
     double pivotValue;
     
     pivotValue = gcv_array[pivotIndex];

     /* Move pivot to end */
     swap_scalar(&(gcv_array[pivotIndex]), &(gcv_array[right]));
     swap_sp(&(par_array[pivotIndex]), &(par_array[right]));

     storeIndex = left - 1;

     for(i = left; i <= right-1; i++)
         if (gcv_array[i] <= pivotValue)
         {
             storeIndex++;
             swap_scalar(&(gcv_array[storeIndex]), &(gcv_array[i]));
             swap_sp(&(par_array[storeIndex]), &(par_array[i]));
         }

     /* Move pivot to its final place */
     swap_scalar(&(gcv_array[right]), &(gcv_array[storeIndex + 1]));
     swap_sp(&(par_array[right]), &(par_array[storeIndex + 1]));

     return (storeIndex + 1);
}

void qsort_knot_psgen(double * gcv_array, struct SP_1D ** par_array, int left, int right, int numSort)
{
     /* Only sorts first numSort items */

     int pivotIndex, pivotNewIndex;

     if (right > left)
     {
           pivotIndex = left;
           pivotNewIndex = qsort_knot_partition_psgen(gcv_array, par_array, 
                         left, right, pivotIndex);
           qsort_knot_psgen(gcv_array, par_array, left, pivotNewIndex - 1, numSort);
           if (pivotNewIndex < numSort)
                 qsort_knot_psgen(gcv_array, par_array, pivotNewIndex + 1, right, numSort);
     }
}
