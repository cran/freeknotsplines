#ifndef PSPLINE_LIB
#define PSPLINE_LIB

#include <stdlib.h> 
#include <math.h>
#include <string.h>
#include <float.h>

#define CALLOC calloc

struct SP_1D /* 1D Spline Structure */ {
	int			order;
	int         n;
	double *	coef;
	double *    knot;
};

struct L2_1D_DATA{  /* Least Squares Data */
	int			nobs;
	double *    xdata;
	double *    ydata;
};

void interv ( double * xt, int lxt, double x, int * left, int * mflag );

void intervl ( double * xt, int lxt, double x, int * right, int * mflag );

double bvalue ( double *t, double *bcoef, int n, int k, double x, int jderiv );

double bvaluel ( double *t, double *bcoef, int n, int k, double x, int jderiv );

void bsplvb (double * t, int jhigh, int index, double x, 
			 int left, double * biatx );

struct SP_1D * sp_1d_initialize(int order, int n,
                                 int  * inform);

double sp_1d_value(struct SP_1D * sp, double x, int deriv);

double ** get_mat(int n1, int n2);

void sp_1d_set_knots(struct SP_1D * sp, double a, double b);

void free_mat(double ** mat);

double plus(double x, int k);

void  bchslv ( double **w, int nbands, int nrow, double *b );

void bchfac ( double **w , int nbands, int nrow, double *diag );

void  get_L2_1D_spline( struct L2_1D_DATA * data, struct SP_1D * sp );

int min_int(int i, int j);

struct L2_1D_DATA * data_1d_initialize(int nobs, int  * inform) ;

struct SP_1D * spline_1d_copy(struct SP_1D * sp);

void free_L2_1D_DATA(struct L2_1D_DATA *p);

void free_SP_1D(struct SP_1D * s);

void get_random_knots(struct SP_1D * sp, struct L2_1D_DATA * data);

void gnomesort(int n, double * ar);

double get_chg_basis_mat(double* t, int n, int k, int i, int j);

double get_omega(double* t, int n, int k, int i, int j);

void  get_pen_1D_spline( struct L2_1D_DATA * data, struct SP_1D * sp , double lambda);

void bchinvb( double **w , int nbands, int nrow, double *diag );

double trace_hat_matrix( struct L2_1D_DATA * data, struct SP_1D * sp , double lambda);

double trace_hat_matrix_fit( struct L2_1D_DATA * data, struct SP_1D * sp , double lambda);

double get_L2_error(struct SP_1D * sp, struct L2_1D_DATA * data);

double gcv( struct L2_1D_DATA * data, struct SP_1D * sp, double lambda, int minspace);

double gcv_fit(struct L2_1D_DATA * data, struct SP_1D * sp, double lambda, int minspace);

int check_knots(struct L2_1D_DATA * data, struct SP_1D * sp, int minspace);

void swap_scalar(double* x1, double* x2);

void swap_sp(struct SP_1D ** sp1, struct SP_1D ** sp2);

int qsort_data_partition(struct L2_1D_DATA * out, int left, int right, int pivotIndex);

void qsort_data(struct L2_1D_DATA * out, int left, int right);

void sort_data(struct L2_1D_DATA * in, struct L2_1D_DATA * out);

void scale_data(struct L2_1D_DATA * in, struct L2_1D_DATA * out, double * shift, double * scale);

double GSJS(struct L2_1D_DATA * data);

double get_chg_basis_mat_order_2(double* t, int i, int j);

double get_chg_basis_mat_order_3(double* t, int i, int j);

double get_chg_basis_mat_order_4(double* t, int i, int j);

#endif
