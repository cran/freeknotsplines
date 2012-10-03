#include "lecuyer.h"

void RandomInitialise(int stream, int nstreams, int seed);
double RandomUniform(void);
double RandomGaussian(double mean,double stddev);
int RandomInt(int lower,int upper);
double RandomDouble(double lower,double upper);

RngStream currstream;

void RandomInitialise(int stream, int nstreams, int seed)
{ 
    unsigned long seedvec[6], i;

    /* set seed */
    for (i = 0; i <= 5; i++)
         seedvec[i] =  (unsigned long) seed;

    RngStream_SetPackageSeed(seedvec);
    
    /* create stream */
    currstream = RngStream_CreateStream("");
    RngStream_ResetStartStream(currstream);
   
    for (i = 0; i < stream; i++)
         currstream = RngStream_CreateStream("");  
}

double RandomUniform(void)
{
      /* generate double precision random number */
      return RngStream_RandU01(currstream);
}


/*
  ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
  The function returns a normally distributed pseudo-random number
  with a given mean and standard devaiation.  Calls are made to a
  function subprogram which must return independent random
  numbers uniform in the interval (0,1).
  The algorithm uses the ratio of uniforms method of A.J. Kinderman
  and J.F. Monahan augmented with quadratic bounding curves.
*/
double RandomGaussian(double mean,double stddev)
{
   double  q,u,v,x,y;

	/*  
		Generate P = (u,v) uniform in rect. enclosing acceptance region 
      Make sure that any random numbers <= 0 are rejected, since
      gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0.
	*/
   do {
      u = RandomUniform();
      v = RandomUniform();
   	if (u <= 0.0 || v <= 0.0) {
       	u = 1.0;
       	v = 1.0;
   	}
      v = 1.7156 * (v - 0.5);

      /*  Evaluate the quadratic form */
      x = u - 0.449871;
   	y = fabs(v) + 0.386595;
      q = x * x + y * (0.19600 * y - 0.25472 * x);

      /* Accept P if inside inner ellipse */
      if (q < 0.27597)
			break;

      /*  Reject P if outside outer ellipse, or outside acceptance region */
    } while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u));

    /*  Return ratio of P's coordinates as the normal deviate */
    return (mean + stddev * v / u);
}

/*
   Return random integer within a range, lower -> upper INCLUSIVE
*/
int RandomInt(int lower,int upper)
{
   return((int)(RandomUniform() * (upper - lower + 1)) + lower);
}

/*
   Return random float within a range, lower -> upper
*/
double RandomDouble(double lower,double upper)
{
   return((upper - lower) * RandomUniform() + lower);
}

