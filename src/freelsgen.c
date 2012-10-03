#include "pspline_lib.h"
#include "lecuyer.h"

struct L2_1D_DATA *  dat_lsgen;
struct SP_1D * sp_lsgen;

const int POP_SIZE_LSGEN = 1000;
const int NUMBER_FIT_LSGEN = 100;
const int NUMBER_MUTATE_LSGEN = 1000;
const int NUMBER_ITER_LSGEN = 200;
const int MIN_SPACE_LSGEN = 2;
const int NUMBER_LOOPS_LSGEN = 20;

const double END_EPSILON_LSGEN = 0.001;
const int END_CHECK_LSGEN = 10;

void mutate_lsgen(struct SP_1D * child);
void reproduce_lsgen(struct SP_1D * par1, struct SP_1D * par2, struct SP_1D * child);

double knot_loop_lsgen(double shift, double scale);

int qsort_knot_partition_lsgen(double * gcv_array, struct SP_1D ** par_array, 
         int left, int right, int pivotIndex);
void qsort_knot_lsgen(double * gcv_array, struct SP_1D ** par_array, 
         int left, int right, int numSort);

void slave_lsgen(double shift, double scale, int seed, int stream);

void freelsgen(int* sampsize, double* xdata, double* ydata, int* ord,
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
dat_lsgen = data_1d_initialize(nobs, &inform);
sort_data(origdat, sortdat);
scale_data(sortdat, dat_lsgen, &shift, &scale);

sp_lsgen = sp_1d_initialize( order, n,  &inform);
sp_1d_set_knots(sp_lsgen, 0., 1.);

best_sp = sp_1d_initialize( order, n,  &inform);
sp_1d_set_knots(best_sp, 0., 1.);

slave_lsgen(shift, scale, seed, stream);

/* Save information */
/* Optimal knots */ 
for (i = 0; i <*numknot; i++)
    optknot[i] = shift + scale * sp_lsgen->knot[i+*ord];

/* trace of hat matrix */
*trace_hat = trace_hat_matrix(dat_lsgen, sp_lsgen, 0.0);

/* GCV */
*GCV = gcv(dat_lsgen, sp_lsgen, 0.0, MIN_SPACE_LSGEN);

/* GSJS */
*GSJS_value = GSJS(dat_lsgen);
}

void slave_lsgen(double shift, double scale, int seed, int stream)
{
     int i;
     double currgcv, bestgcv;

     struct SP_1D  * best_sp;       
     best_sp = spline_1d_copy(sp_lsgen);
     
     RandomInitialise(stream, stream + 1, seed);
     sp_1d_set_knots(sp_lsgen, 0., 1.);
     sp_1d_set_knots(best_sp, 0., 1.);
     bestgcv = INFINITY;

     for (i = 0; i < NUMBER_LOOPS_LSGEN; i++)
     {
          currgcv = knot_loop_lsgen(shift, scale);
          if (currgcv < bestgcv)
          {
              free_SP_1D(best_sp);
              best_sp = spline_1d_copy(sp_lsgen);
              bestgcv = currgcv;
          }
     }

     free_SP_1D(sp_lsgen);
     sp_lsgen = spline_1d_copy(best_sp);
     get_L2_1D_spline(dat_lsgen, sp_lsgen);
     free_SP_1D(best_sp);
}

double knot_loop_lsgen(double shift, double scale)
{
int i,j;
int inform;
int rand1, rand2;

struct SP_1D  * best_sp;
struct SP_1D  * par_array[POP_SIZE_LSGEN];
struct SP_1D  * child_array[POP_SIZE_LSGEN];

double  gcv_array[POP_SIZE_LSGEN];

double  bestgcv; 

double old_crit, new_crit, conv_crit;
int end_flag = 0;

best_sp = spline_1d_copy(sp_lsgen);
bestgcv = INFINITY;

old_crit = INFINITY;

/* initialize parent splines */
for (i = 0; i < POP_SIZE_LSGEN; i++)
    par_array[i] = sp_1d_initialize( sp_lsgen->order, sp_lsgen->n,  &inform);

for (i = 0; i < POP_SIZE_LSGEN; i++)
    sp_1d_set_knots(par_array[i], 0., 1.);

/* generate random splines for parents */
for (i = 0; i < POP_SIZE_LSGEN; i++)
    get_random_knots(par_array[i], dat_lsgen);

/* Begin major loop  */

for(j=0; j< NUMBER_ITER_LSGEN; j++){

/* solve L2 problem for parents and compute errors for parents */
for (i = 0; i < POP_SIZE_LSGEN; i++)
{
    get_L2_1D_spline( dat_lsgen, par_array[i]);
    gcv_array[i] = gcv( dat_lsgen, par_array[i], 0.0, MIN_SPACE_LSGEN);
}

/* sort parent splines by psse */
/* note: only the fittest splines need to be sorted */
qsort_knot_lsgen(gcv_array, par_array, 0, POP_SIZE_LSGEN-1, NUMBER_FIT_LSGEN);

/* update best spline */
if(gcv_array[0] < bestgcv){
	free_SP_1D(best_sp);
	bestgcv = gcv_array[0];
	best_sp = spline_1d_copy(par_array[0]);
}

/* update termination criterion */
if (j % END_CHECK_LSGEN == 0)
{
     new_crit = gcv_array[0] - get_L2_error(best_sp, dat_lsgen)/(dat_lsgen->nobs);
     conv_crit  = fabs(new_crit - old_crit)/old_crit;
     if (conv_crit < END_EPSILON_LSGEN)
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
for (i = 0; i < POP_SIZE_LSGEN; i++)
     child_array[i] = sp_1d_initialize( sp_lsgen->order, sp_lsgen->n,  &inform);

for (i = 0; i < POP_SIZE_LSGEN; i++)
     sp_1d_set_knots(child_array[i], 0., 1.);

/* create child splines */
for (i = 0; i < POP_SIZE_LSGEN; i++)
{
     rand1 = RandomInt(0, NUMBER_FIT_LSGEN - 1);
     rand2 = RandomInt(0, NUMBER_FIT_LSGEN - 1);
     reproduce_lsgen(par_array[rand1], par_array[rand2], child_array[i]);
}

/* mutate_lsgen child splines */
for (i = 0; i < NUMBER_MUTATE_LSGEN; i++)
{
     rand1 = RandomInt(0, POP_SIZE_LSGEN - 1);
     mutate_lsgen(child_array[rand1]);
}

/* free parent splines */
for (i = 0; i < POP_SIZE_LSGEN; i++)
    free_SP_1D(par_array[i]);

/* copy child splines to parent splines */
for (i = 0; i < POP_SIZE_LSGEN; i++)
    par_array[i] = spline_1d_copy(child_array[i]);

/* free child splines */
for (i = 0; i < POP_SIZE_LSGEN; i++)
    free_SP_1D(child_array[i]);

} /* j loop completed*/

/* free parent splines */
for (i = 0; i < POP_SIZE_LSGEN; i++)
    free_SP_1D(par_array[i]);

free_SP_1D(sp_lsgen);
sp_lsgen = spline_1d_copy(best_sp);

return bestgcv;
}

void reproduce_lsgen(struct SP_1D * par1, struct SP_1D * par2, struct SP_1D * child)
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

void mutate_lsgen(struct SP_1D * child)
{
     int knot_number;

     knot_number = RandomInt(child->order, child->n - 1);
     child->knot[knot_number] = RandomDouble(child->knot[knot_number-1], child->knot[knot_number+1]);
}

int qsort_knot_partition_lsgen(double * gcv_array, struct SP_1D ** par_array, int left, int right, int pivotIndex)
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

void qsort_knot_lsgen(double * gcv_array, struct SP_1D ** par_array, int left, int right, int numSort)
{
     /* Only sorts first numSort items */

     int pivotIndex, pivotNewIndex;

     if (right > left)
     {
           pivotIndex = left;
           pivotNewIndex = qsort_knot_partition_lsgen(gcv_array, par_array, 
                         left, right, pivotIndex);
           qsort_knot_lsgen(gcv_array, par_array, left, pivotNewIndex - 1, numSort);
           if (pivotNewIndex < numSort)
                 qsort_knot_lsgen(gcv_array, par_array, pivotNewIndex + 1, right, numSort);
     }
}
