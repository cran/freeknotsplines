#include "pspline_lib.h"

int ilor = 1;
int ilol = 1;

double RandomUniform(void);

double bvalue ( double *t, double *bcoef, int n, int k, double x, int jderiv )
{
/*  from  * a practical guide to splines *  by c. de boor    
calls  interv
c
calculates value at  x  of  jderiv-th derivative of spline from b-repr.
c  the spline is taken to be continuous from the right, EXCEPT at the
c  rightmost knot, where it is taken to be continuous from the left.
c
c******  i n p u t ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
c******  m e t h o d  ******
c     The nontrivial knot interval  (t[i),t[i+1])  containing  x  is lo-
c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef[.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef[.),                     ,  j .eq. 0
c                   /
c    bcoef[.,j)  =  / bcoef[.,j-1] - bcoef[.-1,j-1]
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t[.+k-j] - t[.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t[.))*a(.) + (t[.+m-1] - x)*a(.-1]
c    a(.,x)  =    ---------------------------------------
c                 (x - t[.))      + (t[.+m-1] - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c 
*/
      const int kmax = 20;
      int   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmj,km1,mflag
                    ,nmi,jdrvp1;

      double   aj[21],dl[21],dr[21], fkmj, bvalue_ret;

      bvalue_ret = 0.;
      if (jderiv >= k)                goto label_99;
/*
c
c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t[i) .lt. t[i+1]   and
c      t[i) .le. x .lt. t[i+1] . If no such i can be found,  x  lies
c      outside the support of  the spline  f , hence  bvalue = 0.
c      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c      at  t[n+k] where it is leftcontinuous.) */
      interv ( t, n+k, x, &i, &mflag );
      if (mflag != 0)                 goto label_99;
/*
c  *** if k = 1 (and jderiv = 0), bvalue = bcoef[i].
*/
      km1 = k - 1;
      if (km1 > 0)                   goto label_1;
      bvalue_ret = bcoef[i];
                                        goto label_99;
/*
c
c  *** store the k b-spline coefficients relevant for the knot interval
c     (t[i),t[i+1]) in aj[1],...,aj[k) and compute dl[j] = x - t[i+1-j],
c     dr[j] = t[i+j] - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t[1] or
c     to t[n+k) appropriately.
*/
label_1:
	  jcmin = 1;
      imk = i - k;
      if (imk >= 0)                   goto label_8;
      jcmin = 1 - imk;
	  for(j=1; j<=i; j++)
		  dl[j] = x - t[i+1-j];
      for(j=i; j<= km1; j++)
          aj[k-j] = 0., dl[j] = dl[i];
    
      
                                        goto label_10;
label_8:
	  for(j=1; j<=km1; j++)
		dl[j] = x - t[i+1-j];
	
label_10: 
	  jcmax = k;
      nmi = n - i;
      if (nmi >= 0)                   goto label_18;
      jcmax = k + nmi;
	  for(j=1; j<= jcmax; j++)
		  dr[j] = t[i+j] - x;

      for(j=jcmax; j<=km1; j++)
		  aj[j+1] = 0.,
	      dr[j] = dr[jcmax];

      
                                        goto label_20;
   
label_18:	  
		for(j=1; j<= km1; j++)
		  dr[j] = t[i+j] - x;								
										
   
label_20:    
		for(jc=jcmin; jc<=jcmax; jc++)
           aj[jc] = bcoef[imk + jc];
		  
		   
/*		   
c
c               *** difference the coefficients  jderiv  times.
*/
      if (jderiv == 0)                goto label_30;
	  for(j=1; j<=jderiv; j++){
		  kmj = k-j;
		  fkmj = (double) kmj;
		  ilo = kmj;
		  for(jj=1; jj<=kmj; jj++){
			  aj[jj] = ((aj[jj+1] - aj[jj])/(dl[ilo] + dr[jj]))*fkmj;
              ilo = ilo - 1;
		  }
	  }
/*
c
c  *** compute value at  x  in (t[i),t[i+1]) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj[1],...,aj[k-jderiv).
*/
label_30:  
    if (jderiv == km1)              goto label_39;
      jdrvp1 = jderiv + 1;
    for(j=jdrvp1; j<=km1; j++){
		kmj = k-j;
		ilo = kmj;
		for(jj=1; jj<=kmj; jj++){
			aj[jj] = (aj[jj+1]*dl[ilo] + aj[jj]*dr[jj])/(dl[ilo]+dr[jj]);
			ilo = ilo - 1;
		}
	}
label_39:
	bvalue_ret = aj[1];
label_99:
	return bvalue_ret;
}

double bvaluel ( double *t, double *bcoef, int n, int k, double x, int jderiv )
{
/*  from  * a practical guide to splines *  by c. de boor    
calls  interv
c
calculates value at  x  of  jderiv-th derivative of spline from b-repr.
c  the spline is taken to be continuous from the left, EXCEPT at the
c  leftmost knot, where it is taken to be continuous from the right.
c
c******  i n p u t ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
c******  m e t h o d  ******
c     The nontrivial knot interval  (t[i),t[i+1])  containing  x  is lo-
c  cated with the aid of  intervl . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef[.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef[.),                     ,  j .eq. 0
c                   /
c    bcoef[.,j)  =  / bcoef[.,j-1] - bcoef[.-1,j-1]
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t[.+k-j] - t[.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t[.))*a(.) + (t[.+m-1] - x)*a(.-1]
c    a(.,x)  =    ---------------------------------------
c                 (x - t[.))      + (t[.+m-1] - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c 
*/
      const int kmax = 20;
      int   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmj,km1,mflag
                    ,nmi,jdrvp1;

      double   aj[21],dl[21],dr[21], fkmj, bvalue_ret;

      bvalue_ret = 0.;
      if (jderiv >= k)                goto label_99;
/*
c
c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t[i) .lt. t[i+1]   and
c      t[i) .lt. x .le. t[i+1] . If no such i can be found,  x  lies
c      outside the support of  the spline  f , hence  bvalue = 0.
c      (The asymmetry in this choice of  i  makes  f  leftcontinuous, except
c      at  t[1] where it is rightcontinuous.) */
      intervl ( t, n+k, x, &i, &mflag );

      if (mflag != 0)                 goto label_99;
/*
c  *** if k = 1 (and jderiv = 0), bvalue = bcoef[i].
*/
      km1 = k - 1;
      if (km1 > 0)                   goto label_1;
      bvalue_ret = bcoef[i];
                                        goto label_99;
/*
c
c  *** store the k b-spline coefficients relevant for the knot interval
c     (t[i),t[i+1]) in aj[1],...,aj[k) and compute dl[j] = x - t[i+1-j],
c     dr[j] = t[i+j] - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t[1] or
c     to t[n+k) appropriately.
*/
label_1:
	  jcmin = 1;
      imk = i - k;
      if (imk >= 0)                   goto label_8;
      jcmin = 1 - imk;
	  for(j=1; j<=i; j++)
		  dl[j] = x - t[i+1-j];
      for(j=i; j<= km1; j++)
          aj[k-j] = 0., dl[j] = dl[i];
    
      
                                        goto label_10;
label_8:
	  for(j=1; j<=km1; j++)
		dl[j] = x - t[i+1-j];
	
label_10: 
	  jcmax = k;
      nmi = n - i;
      if (nmi >= 0)                   goto label_18;
      jcmax = k + nmi;
	  for(j=1; j<= jcmax; j++)
		  dr[j] = t[i+j] - x;

      for(j=jcmax; j<=km1; j++)
		  aj[j+1] = 0.,
	      dr[j] = dr[jcmax];

      
                                        goto label_20;
   
label_18:	  
		for(j=1; j<= km1; j++)
		  dr[j] = t[i+j] - x;								
										
   
label_20:    
		for(jc=jcmin; jc<=jcmax; jc++)
           aj[jc] = bcoef[imk + jc];
		  
		   
/*		   
c
c               *** difference the coefficients  jderiv  times.
*/
      if (jderiv == 0)                goto label_30;
	  for(j=1; j<=jderiv; j++){
		  kmj = k-j;
		  fkmj = (double) kmj;
		  ilo = kmj;
		  for(jj=1; jj<=kmj; jj++){
			  aj[jj] = ((aj[jj+1] - aj[jj])/(dl[ilo] + dr[jj]))*fkmj;
              ilo = ilo - 1;
		  }
	  }
/*
c
c  *** compute value at  x  in (t[i),t[i+1]) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj[1],...,aj[k-jderiv).
*/
label_30:  
    if (jderiv == km1)              goto label_39;
      jdrvp1 = jderiv + 1;
    for(j=jdrvp1; j<=km1; j++){
		kmj = k-j;
		ilo = kmj;
		for(jj=1; jj<=kmj; jj++){
			aj[jj] = (aj[jj+1]*dl[ilo] + aj[jj]*dr[jj])/(dl[ilo]+dr[jj]);
			ilo = ilo - 1;
		}
	}
label_39:
	bvalue_ret = aj[1];
label_99:
	return bvalue_ret;
}



void interv ( double * xt, int lxt, double x, int * left, int * mflag )
{
/*
c  from  * a practical guide to splines *  by C. de Boor    
computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
c
c******  i n p u t  ******
c  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the right, but, by returning  mflag = 0  even if
C        x = xt(lxt), there is the option of having the computed pp function
c        continuous from the left at  xt(lxt) .
c
c******  m e t h o d  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  global varia-
c  ble  ilor . A first check ascertains that  ilor .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilor) .le. x .lt. xt(ilor+1), we set  left =
c  ilor  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilor
c  while also moving  ilor  and  ihi  in the direction of  x , until
c                      xt(ilor) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilor+1 = ihi .
c  left = ilor  is then returned.
c
*/
      int  ihi, istep, middle;
         
      if(ilor <1 || ilor>lxt) ilor = 1;
      ihi = ilor + 1;
      if (ihi < lxt)                  goto label_20;
         if (x >= xt[lxt])            goto label_110;
         if (lxt <= 1)                goto label_90;
         ilor = lxt - 1;
         ihi = lxt;

label_20:
      if (x >= xt[ihi])               goto label_40;
      if (x >= xt[ilor])               goto label_100;
/*
              **** now x .lt. xt(ilor) . decrease  ilor  to capture  x .
*/
      istep = 1;
label_31:
         ihi = ilor;
         ilor = ihi - istep;
         if (ilor <= 1)                  goto label_35;
         if (x >= xt[ilor])              goto label_50;
         istep = istep*2;
										goto label_31;
label_35:
	  ilor = 1;
      if (x < xt[1])                    goto label_90;
                                        goto label_50;
/*
              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
*/
label_40:   
   istep = 1;
label_41:
         ilor = ihi;
         ihi = ilor + istep;
         if (ihi >= lxt)              goto label_45;
         if (x < xt[ihi])             goto label_50;
         istep = istep*2;
                                      goto label_41;
label_45:
      if (x >= xt[lxt])               goto label_110;
      ihi = lxt;
/*
c           **** now xt(ilor) .le. x .lt. xt(ihi) . narrow the interval.
*/
label_50:
      middle = (ilor + ihi)/2;
      if (middle == ilor)              goto label_100;
/*
c     note. it is assumed that middle = ilor in case ihi = ilor+1 .
*/
      if (x < xt[middle])              goto label_53;
         ilor = middle;
                                       goto label_50;
label_53:
      ihi = middle;
                                       goto label_50;
/* c**** set output and return.  */
label_90:
      *mflag = -1;
      *left = 1;
                                        return;
label_100:  
	  *mflag = 0;
      *left = ilor;
                                        return;
label_110:  
	  *mflag = 1;
      if (x == xt[lxt]) *mflag = 0;
      *left = lxt;
label_111:
      if (*left == 1)                  return;
      *left = *left - 1;
      if (xt[*left] < xt[lxt])        return;
          goto label_111;
}


void intervl ( double * xt, int lxt, double x, int * left, int * mflag )
{
/*
c******  i n p u t  ******
c  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   1      0      if   xt(1)  .eq. x 
c   i      0      if   xt(i)  .lt. x .le. xt(i+1)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the left, but, by returning  mflag = 0  even if
C        x = xt(lxt), there is the option of having the computed pp function
c        continuous from the rightt at  xt(1) .
c
c******  m e t h o d  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  global  varia-
c  ble  ilol . A first check ascertains that  ilol .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilol) .lt. x .le. xt(ilol+1), we set  left =
c  ilol  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilol
c  while also moving  ilol  and  ihi  in the direction of  x , until
c                      xt(ilol) .lt. x .le. xt(ihi) ,
c  after which we use bisection to get, in addition, ilol+1 = ihi .
c  left = ilol  is then returned.
c
*/
      int  ihi, istep, middle;
      
      if(ilol <1 || ilol>lxt) ilol = 1;
      ihi = ilol + 1;
      if (ihi < lxt)                  goto label_20;
         if (x >= xt[lxt])            goto label_110;
         if (lxt <= 1)                goto label_90;
         ilol = lxt - 1;
         ihi = lxt;

label_20:
      if (x > xt[ihi])               goto label_40;
      if (x > xt[ilol])               goto label_100;
/*
              **** now x .le. xt(ilol) . decrease  ilol  to capture  x .
*/
      istep = 1;
label_31:
         ihi = ilol;
         ilol = ihi - istep;
         if (ilol <= 1)                  goto label_35;
         if (x > xt[ilol])              goto label_50;
         istep = istep*2;
										goto label_31;
label_35:
	  ilol = 1;
      if (x < xt[1])                    goto label_90;
                                        goto label_50;
/*
              **** now x .gt. xt(ihi) . increase  ihi  to capture  x .
*/
label_40:   
   istep = 1;
label_41:
         ilol = ihi;
         ihi = ilol + istep;
         if (ihi >= lxt)              goto label_45;
         if (x <= xt[ihi])             goto label_50;
         istep = istep*2;
                                      goto label_41;
label_45:
      if (x >= xt[lxt])               goto label_110;
      ihi = lxt;
/*
c           **** now xt(ilol) .lt. x .le. xt(ihi) . narrow the interval.
*/
label_50:
      middle = (ilol + ihi)/2;
      if (middle == ilol)              goto label_100;
/*
c     note. it is assumed that middle = ilol in case ihi = ilol+1 .
*/
      if (x <= xt[middle])              goto label_53;
         ilol = middle;
                                       goto label_50;
label_53:
      ihi = middle;
                                       goto label_50;
/* c**** set output and return.  */
label_90:
      *mflag = -1;
      *left = 1;
                                        return;
label_100:  
	  *mflag = 0;
      *left = ilol;
                                        return;
label_110:  
	  *mflag = 1;
      if (x == xt[lxt]) *mflag = 0;
      *left = lxt;
label_111:
      if (*left == 1)                  return;
      *left = *left - 1;
      if (xt[*left] < xt[lxt])        return;
          goto label_111;
}


void bsplvb (double * t, int jhigh, int index, double x, int left, double * biatx )
{
/*
c  from  * a practical guide to splines *  by c. de boor    
calculates the value of all possibly nonzero b-splines at  x  of order
c
c               jout  =  max( jhigh , (j+1)*(index-1) )
c
c  with knot sequence  t .
c
c******  i n p u t  ******
c  t.....knot sequence, of length  left + jout  , assumed to be nonde-
c        creasing.  a s s u m p t i o n . . . .
c                       t(left)  .lt.  t(left + 1)   .
c   d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
c  jhigh,
c  index.....integers which determine the order  jout = max(jhigh,
c        (j+1)*(index-1))  of the b-splines whose values at  x  are to
c        be returned.  index  is used to avoid recalculations when seve-
c        ral columns of the triangular array of b-spline values are nee-
c        ded (e.g., in  bsplpp  or in  bsplvd ). precisely,
c                     if  index = 1 ,
c        the calculation starts from scratch and the entire triangular
c        array of b-spline values of orders 1,2,...,jhigh  is generated
c        order by order , i.e., column by column .
c                     if  index = 2 ,
c        only the b-spline values of order  j+1, j+2, ..., jout  are ge-
c        nerated, the assumption being that  biatx , j , deltal , deltar
c        are, on entry, as they were on exit at the previous call.
c           in particular, if  jhigh = 0, then  jout = j+1, i.e., just
c        the next column of b-spline values is generated.
c
c  w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
c        posed arbitrarily by the dimension statement for  deltal  and
c        deltar  below, but is  n o w h e r e  c h e c k e d  for .
c
c  x.....the point at which the b-splines are to be evaluated.
c  left.....an integer chosen (usually) so that
c                  t(left) .le. x .le. t(left+1)  .
c
c******  o u t p u t  ******
c  biatx.....array of length  jout , with  biatx(i)  containing the val-
c        ue at  x  of the polynomial of order  jout  which agrees with
c        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
c        t(left+1)) .
c
c******  m e t h o d  ******
c  the recurrence relation
c
c                       x - t(i)              t(i+j+1) - x
c     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
c                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
c
c  is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
c  ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
c  b(left,j)(x), storing the new values in  biatx  over the old. the
c  facts that
c            b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)
c  and that
c            b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)
c  are used. the particular organization of the calculations follows al-
c  gorithm  (8)  in chapter x of the text.
c
      
current fortran standard makes it impossible to specify the length of
c  t  and of  biatx  precisely without the introduction of otherwise
c  superfluous additional arguments.
*/


/* NOTE - THE FOLLOWING VARIABLES (j, deltal, and deltar) WERE ORIGINALLY DECLARED STATIC. 
   FOR ALL THE SPLINE-FITTING ALGORITHMS, THIS CODE IS CALLED WITH INDEX = 1.  THEREFORE, 
   THERE IS NO NEED FOR THEM TO BE DECLARED STATIC.   FURTHERMORE, SINCE OPENMP IS A 
    SHARED MEMORY SYSTEM, STATIC VARIABLES MAY BE SHARED BETWEEN THREADS.  THEREFORE, 
   THEY ARE NO LONGER STATIC VARIABLES. */ 

      int  j = 0;
	 double deltal[21],deltar[21]; 

	int jmax = 20;
      int   i, jp1;
      double  saved, term;          
      
	  if(j <= 0) j=1;	
      if(index == 1) 
		  goto label_10;
	  else
		  goto label_20;

label_10:                          
      j = 1;
      biatx[1] = 1.;
      if (j >= jhigh)                 goto label_99;

label_20:
         jp1 = j + 1;
         deltar[j] = t[left+j] - x;
         deltal[j] = x - t[left+1-j];
         saved = 0.;
		 for(i=1; i<=j; i++){
            term = biatx[i]/(deltar[i] + deltal[jp1-i]);
            biatx[i] = saved + deltar[i]*term;
            saved = deltal[jp1-i]*term;
		 }
         biatx[jp1] = saved;
         j = jp1;
         if (j < jhigh)              goto label_20;
label_99:
      return;
}
    
struct SP_1D * sp_1d_initialize(int order, int n,
                                 int  * inform)
{
/*  Purpose: This function allocates and initializes all the pointers
    in the SP_1D data structure.  

    Return Value:  A pointer to the newly initialized SP_1D structure.


    Arguments:

        order           Order of spline (input)
        n               Spline space dimension (input)
        
        
        inform          Error code (output)


    Fatal Errors (which return a NULL pointer):

        inform          Error Condition

          1             order  <= 0
          2             n      <= 0
          4             allocation failure for some member of SP_2D_TP

    */

    struct SP_1D * sp;
    double * temp_ptr1;
    

    /* Initial error testing */
    *inform = 0;
    if(order <= 0)
        *inform = 1;
    else if(n <= 0)
        *inform = 2;
    if(*inform > 0) return NULL;



            /*  Allocate space for the SP_1D structure  */
    sp = (struct SP_1D *) CALLOC(1, sizeof(struct SP_1D));
    if((void *)sp == NULL){
        *inform = 4;
        return NULL;
    }


            /*  Place numerical values in SP_1D structure  */
    sp->order       = order;
    sp->n		    = n;  

            /*  Allocate space for the coefficients
                and set the SP_1D pointers  */
    temp_ptr1   = (double *) CALLOC(n, sizeof(double));
    if((void *)temp_ptr1 == NULL){
        *inform = 4;
        return NULL;
    }

    
    
    (sp->coef)= temp_ptr1;

            /*  Allocate space for the knots
                and set the SP_1D pointers  */
   
    temp_ptr1   = (double *) CALLOC(n+order , sizeof(double));
    if((void *)temp_ptr1 == NULL){
        *inform = 4;
        return NULL;}

	sp->knot = temp_ptr1;

    
	
	return sp;

}

void sp_1d_set_knots(struct SP_1D * sp, double a, double b)
{
/*  Purpose: This function returns equally spaced knots for the 1d 
             spline.  The x-knots partition the interval [a, b].

    Return Value:  void


    Arguments:

        sp              A 1d spline (input)
        a               Left end point for the x-knots (input)
        b               Right end point for the x-knots (input)
   
      
    There are no Fatal Errors 

    */
int i;
double temp;

for(i=0; i<sp->order; i++){
	sp->knot[i] = a;
    sp->knot[sp->n+i] = b;
}

temp = b-a;
for(i=sp->order; i<sp->n; i++)
	sp->knot[i] = a + ((i-sp->order+1)*temp)/(sp->n - sp->order+1);
}


double sp_1d_value(struct SP_1D * sp, double x, int deriv)
{ 
/*  Purpose: This function computes the values of a 1d spline
             or one of its derivatives at x. 

    Return Value:  The value of the (deriv)-th derivative of the
	               the spline stored in sp at x. 



    Arguments:

        sp				1d spline (input)
		x               vector of x-values(input)
        y               vector of y-values (input)
		deriv			the order of the derivative of sp at x (input)
		                deriv must be non-negative.
        
        
        There are no Fatal Errors.

        

    */	
	return bvalue(sp->knot -1, sp->coef-1, sp->n, sp->order,
		          x, deriv);
}
	
double ** get_mat(int n1, int n2)
{
/*  Purpose: This function allocates space for a n1 by n2 doubly
             indexed array of doubles.

    Return Value:  n1 by n2 matrix implemented as pointer to pointer 
	               to double.  The indexing starts at [0][0] to 
				   [n1-1][n2-1].  And the storage is allocated contiguously.
	               
				   

    Arguments:

        n1				row index (input)
		n2              column index (input)
        
        
        There are no Fatal Errors.

	Related Functions:  free_mat

 */
	int i;
	double ** mat, * temp_ptr;

            /*  Allocate space for the array  */
    temp_ptr   = (double *) CALLOC(n1*n2, sizeof(double));
    if((void *)temp_ptr == NULL){
       /* *inform = 4; */
        return NULL;
    }

	mat = (double **) CALLOC(n1, sizeof(double *));
    if((void *)(mat) == NULL){
        /**inform = 4;*/
        return NULL;
    }

    for(i=0; i< n1; i++)
        mat[i] = &(temp_ptr[i*n2]);

	return mat;
}

void free_mat(double ** mat)
{
/*  Purpose: This function frees space allocated by get_mat.

    Return Value:  void
	               
				   

    Arguments:

        mat				a matrix allocated by get_mat (input)
        
        
        There are no Fatal Errors.

	Related Functions:  get_mat

 */
	free(*mat);
	free(mat);
}

double plus(double x, int k)
{
/*  Purpose: This function returns (x)_+^k,  that is, 
             x^k if x>0 or 0 if x<= 0; 


    Return Value:  (x)_+^k
				   
    Arguments:

        x				argument for the "+" function (input)
		k               power for the "+" function (input)
        
    
        There are no Fatal Errors.

 */
	if(x > 0) return pow(x, k);
	return 0.;
}


void  bchslv ( double **w, int nbands, int nrow, double *b )
{/*
c  from  * a practical guide to splines *  by c. de boor    
c  solves the linear system     c*x = b   of order  n r o w  for  x
c  provided  w  contains the cholesky factorization for the banded (sym-
c  metric) positive definite matrix  c  as constructed in the subroutine
c    b c h f a c  (quo vide).
c
c******  i n p u t  ******
c  nrow.....is the order of the matrix  c .
c  nbands.....indicates the bandwidth of  c .
c  w.....contains the cholesky factorization for  c , as output from
c        subroutine bchfac  (quo vide).
c  b.....the vector of length  n r o w  containing the right side.
c
c******  o u t p u t  ******
c  b.....the vector of length  n r o w  containing the solution.
c
c******  m e t h o d  ******
c  with the factorization  c = l*d*l-transpose  available, where  l  is
c  unit lower triangular and  d  is diagonal, the triangular system
c  l*y = b  is solved for  y (forward substitution), y is stored in  b,
c  the vector  d**(-1)*y is computed and stored in  b, then the triang-
c  ular system  l-transpose*x = d**(-1)*y is solved for  x (backsubstit-
c  ution).
*/
      int  j, jmax, n, nbndm1;
      if (nrow > 1)                  goto label_21;
      b[1] = b[1]*w[1][1];
                                        return;

/*
c                 forward substitution. solve l*y = b for y, store in b.
*/
label_21: nbndm1 = nbands - 1;
     
	for(n=1; n<=nrow; n++){
         jmax = min_int(nbndm1,nrow-n);
         if (jmax < 1)               goto label_30;
		 for(j=1; j<= jmax; j++){
            b[j+n] = b[j+n] - w[j+1][n]*b[n];
		}
label_30:     n=n;
		  }
/*
c     backsubstitution. solve l-transp.x = d**(-1)*y  for x, store in b.
*/
			n = nrow;    
label_39:   b[n] = b[n]*w[1][n];  
         jmax = min_int(nbndm1,nrow-n);
         if (jmax < 1)               goto label_40;
			 for(j=1; j<=jmax; j++){
			      b[n] = b[n] - w[j+1][n]*b[j+n];
			 }
label_40:    n = n-1; 
         if (n > 0)                  goto label_39;
return;
}


void bchfac ( double **w , int nbands, int nrow, double *diag )
{
/*
c  from  * a practical guide to splines *  by c. de boor    
constructs cholesky factorization
c                     c  =  l * d * l-transpose
c  with l unit lower triangular and d diagonal, for given matrix c of
c  order  n r o w , in case  c  is (symmetric) positive semidefinite
c  and  b a n d e d , having  n b a n d s  diagonals at and below the
c  main diagonal.
c
c******  i n p u t  ******
c  nrow.....is the order of the matrix  c .
c  nbands.....indicates its bandwidth, i.e.,
c          c(i,j) = 0 for i-j .ge. nbands .
c  w.....workarray of size (nbands,nrow)  containing the  nbands  diago-
c        nals in its rows, with the main diagonal in row  1 . precisely,
c        w(i,j)  contains  c(i+j-1,j), i=1,...,nbands, j=1,...,nrow.
c          for example, the interesting entries of a seven diagonal sym-
c        metric matrix  c  of order  9  would be stored in  w  as
c
c
c
c
c
c
c        all other entries of  w  not identified in this way with an en-
c        try of  c  are never referenced .
c  diag.....is a work array of length  nrow .
c
c******  o u t p u t  ******
c  w.....contains the cholesky factorization  c = l*d*l-transp, with
c        w(1,i) containing  1/d(i,i)
c        and  w(i,j)  containing  l(i-1+j,j), i=2,...,nbands.
c
c******  m e t h o d  ******
c   gauss elimination, adapted to the symmetry and bandedness of  c , is
c   used .
c     near zero pivots are handled in a special way. the diagonal ele-
c  ment c(n,n) = w[1][n] is saved initially in  diag[n], all n. at the n-
c  th elimination step, the current pivot element, viz.  w[1][n], is com-
c  pared with its original value, diag[n]. if, as the result of prior
c  elimination steps, this element has been reduced by about a word
c  length, (i.e., if w[1][n]+diag[n] .le. diag[n]), then the pivot is de-
c  clared to be zero, and the entire n-th row is declared to be linearly
c  dependent on the preceding rows. this has the effect of producing
c   x(n) = 0  when solving  c*x = b  for  x, regardless of  b. justific-
c  ation for this is as follows. in contemplated applications of this
c  program, the given equations are the normal equations for some least-
c  squares approximation problem, diag[n] = c(n,n) gives the norm-square
c  of the n-th basis function, and, at this point,  w[1][n]  contains the
c  norm-square of the error in the least-squares approximation to the n-
c  th basis function by linear combinations of the first n-1 . having
c  w[1][n]+diag[n] .le. diag[n] signifies that the n-th function is lin-
c  early dependent to machine accuracy on the first n-1 functions, there
c  fore can safely be left out from the basis of approximating functions
c     the solution of a linear system
c                       c*x = b
c   is effected by the succession of the following  t w o  calls:
c     call bchfac ( w, nbands, nrow, diag )       , to get factorization
c     call bchslv ( w, nbands, nrow, b )          , to solve for x.
c
*/
      int  i,imax,j,jmax,n;
      double    ratio;
      if (nrow > 1)                  goto label_9;
      if (w[1][1] > 0.) w[1][1] = 1./w[1][1];
                                        return;
/*      store diagonal of  c  in  diag */
label_9: 
	  
		  for(n=1; n<= nrow; n++){
		    diag[n] = w[1][n];
		  }
/* factorization */
    
	for(n=1; n<= nrow; n++){
         if (w[1][n]+diag[n] > diag[n]) goto label_15;
		for(j=1; j<= nbands; j++){
			w[j][n] = 0.;
		}
                                        goto label_20;
label_15:    w[1][n] = 1./w[1][n];
         imax = min_int(nbands-1,nrow - n);
         if (imax < 1)               goto label_20;
         jmax = imax;
		 for(i=1; i<= imax; i++){
            ratio = w[i+1][n]*w[1][n];
			for(j=1; j<= jmax; j++){
				w[j][n+i] = w[j][n+i] - w[j+i][n]*ratio;
			}
            jmax = jmax - 1;
            w[i+1][n] = ratio;
		 }
label_20:    n=n; 
	}
return;
    
}

void  get_L2_1D_spline( struct L2_1D_DATA * data, struct SP_1D * sp )
{
/* l2appr ( t, n, k, q, diag, bcoef )
c  from  * a practical guide to splines *  by c. de boor    
c  to be called in main program  l 2 m a i n .
calls subprograms  bsplvb, bchfac/slv
c
constructs the (discrete) l2-approximation by splines of order
c  k  with knot sequence  t(1), ..., t(n+k)  to given data points
c  ( tau(i), gtau(i) ), i=1,...,ntau. the b-spline coefficients
c  b c o e f   of the approximating spline are determined from the
c  normal equations using cholesky's method.
c
c******  i n p u t  ******
c  t(1), ..., t(n+k)  the knot sequence
c  n.....the dimension of the space of splines of order k with knots t.
c  k.....the order
c
c  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
c        sed by the arbitrary dimension statement for  biatx  below, but
c        is  n o w h e r e   c h e c k e d   for.
c
c******  w o r k  a r r a y s  ******
c  q....a work array of size (at least) k*n. its first  k  rows are used
c       for the  k  lower diagonals of the gramian matrix  c .
c  diag.....a work array of length  n  used in bchfac .
c
c******  i n p u t  via  c o m m o n  /data/  ******
c  ntau.....number of data points
c  (tau(i),gtau(i)), i=1,...,ntau     are the  ntau  data points to be
c        fitted .
c
c******  o u t p u t  ******
c  bcoef(1), ..., bcoef(n)  the b-spline coeffs. of the l2-appr.
c
c******  m e t h o d  ******
c  the b-spline coefficients of the l2-appr. are determined as the sol-
c  ution of the normal equations
c     sum ( (b(i),b(j))*bcoef(j) : j=1,...,n)  = (b(i),g),
c                                               i = 1, ..., n .
c  here,  b(i)  denotes the i-th b-spline,  g  denotes the function to
c  be approximated, and the  i n n e r   p r o d u c t  of two funct-
c  ions  f  and  g  is given by
c      (f,g)  :=  sum ( f(tau(i))*g(tau(i)) : i=1,...,ntau) .
c  the arrays  t a u  and  w e i g h t  are given in common block
c   d a t a , as is the array  g t a u  containing the sequence
c  g(tau(i)), i=1,...,ntau.
c  the relevant function values of the b-splines  b(i), i=1,...,n, are
c  supplied by the subprogram  b s p l v b .
c     the coeff.matrix  c , with
c           c(i,j)  :=  (b(i), b(j)), i,j=1,...,n,
c  of the normal equations is symmetric and (2*k-1)-banded, therefore
c  can be specified by giving its k bands at or below the diagonal. for
c  i=1,...,n,  we store
c   (b(i),b(j))  =  c(i,j)  in  q(i-j+1,j), j=max(1, i-k+1), ... , i
c  and the right side
c   (b(i), g )  in  bcoef(i) .
c  since b-spline values are most efficiently generated by finding sim-
c  ultaneously the value of  e v e r y  nonzero b-spline at one point,
c  the entries of  c  (i.e., of  q ), are generated by computing, for
c  each ll, all the terms involving  tau(ll)  simultaneously and adding
c  them to all relevant entries.
*/
      
	  const int kmax = 20;
	  int k, n, i, j, jj, left, leftmk, ll, mm, ntau;
	  double * bcoef, *diag, ** q, *t, biatx[21], dw;
	  double *tau, *gtau;
		        
        n = sp->n;
        k = sp->order;
	  t = sp->knot -1;
	  ntau = data->nobs;
	  tau = data->xdata -1;
	  gtau = data->ydata -1;
	  q = get_mat(k+1, n+1);
	  diag =  (double *) CALLOC(n+1, sizeof(double));
        bcoef = (double *) CALLOC(n+1, sizeof(double));
      
        left = k;
        leftmk = 0;
	  for(ll=1; ll<= ntau; ll++){
      
   /*   locate  l e f t  s.t. tau(ll) in (t(left),t(left+1)) */
label_10:   if (left == n)            goto label_15;
            if (tau[ll] < t[left+1])  goto label_15;
            left = left+1;
            leftmk = leftmk + 1;
			goto label_10;
label_15:   bsplvb ( t, k, 1, tau[ll], left, biatx );
/*
c        biatx(mm) contains the value of b(left-k+mm) at tau(ll).
c        hence, with  dw := biatx(mm), the number dw*gtau(ll)
c        is a summand in the inner product
c           (b(left-k+mm), g)  which goes into  bcoef(left-k+mm)
c        and the number biatx(jj)*dw is a summand in the inner product
c           (b(left-k+jj), b(left-k+mm)), into  q(jj-mm+1,left-k+mm)
c        since  (left-k+jj) - (left-k+mm) + 1  =  jj - mm + 1 .
*/
         
		for(mm=1; mm<=k; mm++){
            dw = biatx[mm];
            j = leftmk + mm;
            bcoef[j] = dw*gtau[ll] + bcoef[j];
            i = 1;
         
				for(jj=mm; jj<=k; jj++){
                  q[i][j] = biatx[jj]*dw + q[i][j];
                  i = i + 1;
				}
			 }
	  }
/*
c
c             construct cholesky factorization for  c  in  q , then use
c             it to solve the normal equations
c                    c*x  =  bcoef
c             for  x , and store  x  in  bcoef .
*/
      bchfac ( q, k, n, diag );
      bchslv ( q, k, n, bcoef );
	  for(i=0; i<sp->n; i++) 
		  sp->coef[i] = bcoef[i+1];
	  free_mat(q), free(diag), free(bcoef);
return;
}


int min_int(int i, int j)
{int k;
	k = i;
	if(k > j) k = j;
	return k;
}
	
struct L2_1D_DATA * data_1d_initialize(int nobs, int  * inform) 
                                 
{
/*  Purpose: This function allocates and initializes all the pointers
    in the L2_1D_DATA data structure.  

    Return Value:  A pointer to the newly initialized TP_1D_DATA structure


    Arguments:

        nobs            number of x-values (input)
        
        
        
        inform          Error code (output)


    Fatal Errors (which return a NULL pointer):

        inform          Error Condition

          1             nobs  <= 0
          
          4             allocation failure for some member of L2_1D_DATA

    */

    struct L2_1D_DATA * dat;
	double * temp_ptr1;
    

    /* Initial error testing */
    *inform = 0;
    if(nobs <= 0)
        *inform = 1;
    
    
    if(*inform > 0) return NULL;



            /*  Allocate space for the L2_1D_DATA structure  */
    dat = (struct L2_1D_DATA *) CALLOC(1, sizeof(struct L2_1D_DATA));
    if((void *)dat == NULL){
        *inform = 4;
        return NULL;
    }


            /*  Place numerical values in L2_1D_DATA structure  */
    dat->nobs       = nobs;
        

            /*  Allocate space for the values
                and set the L2_1D_DATA pointers  */
    temp_ptr1   = (double *) CALLOC(2*nobs, sizeof(double));
    if((void *)temp_ptr1 == NULL){
        *inform = 4;
        return NULL;
    }

	dat->xdata  = temp_ptr1; 
	dat->ydata  = temp_ptr1 + nobs;

	
	return dat;

}

struct SP_1D * spline_1d_copy(struct SP_1D * sp)
{
	int i, inform;
	struct SP_1D * copy_sp;

	copy_sp = sp_1d_initialize(sp->order, sp->n, &inform);

	for(i=0; i< sp->order + sp->n ; i++)
		copy_sp->knot[i] = sp->knot[i];
	for(i=0; i< sp->n ; i++)
		copy_sp->coef[i] = sp->coef[i];
	return copy_sp;
}


void free_L2_1D_DATA(struct L2_1D_DATA *p)
{
	/*free(p->ydata); 
	free(p->xdata); 
	
	these statements are unnecessary 
	due to the allocation mechanisim*/
	free(p);
}

void free_SP_1D(struct SP_1D * s)
{
	free(s->coef);
	free(s->knot);
	free(s);
}

void get_random_knots(struct SP_1D * sp, struct L2_1D_DATA * data ){
int i, num_int_knots;
double *  temp, low, high;

num_int_knots = sp->n - sp->order;
temp = (sp->knot + sp->order);

/*  Note this assumes that the xdata field is ordered */
low = data->xdata[0];
high = data->xdata[data->nobs -1];

for(i=0; i<num_int_knots; i++)
	temp[i] = low +RandomUniform()*(high - low);
	
gnomesort(num_int_knots, temp);

/*
for(i=0; i<num_int_knots; i++)
	sp->knot[sp->order + i] = temp[i];

free(temp);
*/
}


void gnomesort(int n, double * ar){
	int i = 0;
	while(i<n){
		if (i==0 || ar[i-1] <= ar[i]) i++;
		else {double tmp = ar[i]; ar[i] = ar[i-1]; ar[--i] = tmp;}
	}

}

double get_chg_basis_mat(double* t, int n, int k, int i, int j)
{
      /* This function returns the coefficient of the jth TPB function in the expansion 
         of the ith B spline.  This is the (j, i)th element of the change of basis matrix
         for B splines in terms of the truncated power basis. */

      double * bcoef;
      double answer;
      double leftderiv, rightderiv;
      int ii;

      /* take advantage of the fact that the change of basis matrix is lower triangular
         and banded */
      if (i > j)
          return 0.0; 
      if (i <= j - (k + 1))
          return 0.0;      


      /* allocate memory */
      bcoef  = (double *) CALLOC(n + 1, sizeof(double));

      for (ii = 0; ii <= n; ii++)
         bcoef[ii] = 0;

      bcoef[i] = 1;

      if (j <= k)
      {
          answer = bvalue(t, bcoef, n, k, t[1], j - 1);
          for (ii = 1; ii <= j - 1; ii++)
               answer /= (double) ii;
      }
      else
      {
          if (k == 2)
              answer = get_chg_basis_mat_order_2(t, i, j);
          else if (k == 3)
              answer = get_chg_basis_mat_order_3(t, i, j);
          else if (k == 4)
              answer = get_chg_basis_mat_order_4(t, i, j);
          else
          {
              leftderiv = bvaluel(t, bcoef, n, k, t[j], k - 1);
              rightderiv = bvalue(t, bcoef, n, k, t[j], k - 1);
              answer = rightderiv - leftderiv;
              for (ii = 1; ii <= k - 1; ii++)
                   answer /= (double) ii;
          }
      }
      
      free(bcoef);
      return answer;
}

double get_omega(double* t, int n, int k, int i, int j)
{

     /* This function returns the (i, j)^th element of the symmetric penalty matrix Omega.
     Here Omega = V^T D_m V where V is the change of basis matrix and D_m is a diagonal matrix 
     with the first m elements equal to 0 and the rest equal to 1. */
 
     double answer = 0.0;
     int ii;
     double entry1, entry2;

     for(ii = k + 1; ii <= n; ii++)
           if ((j <= ii) && (j > ii - (k + 1))) /* otherwise entry2 = 0 */
           {
                 entry1 = get_chg_basis_mat(t, n, k, i, ii);
                 if (entry1 != 0)
                 {
                     entry2 = get_chg_basis_mat(t, n, k, j, ii);
                     answer += entry1 * entry2;
                 }
            }
 
     return(answer);
}

void  get_pen_1D_spline( struct L2_1D_DATA * data, struct SP_1D * sp , double lambda)
{
/* l2appr ( t, n, k, q, diag, bcoef )
c  from  * a practical guide to splines *  by c. de boor    
c  to be called in main program  l 2 m a i n .
calls subprograms  bsplvb, bchfac/slv
c
constructs the (discrete) penalized l2-approximation by splines of order
c  k  with knot sequence  t(1), ..., t(n+k)  to given data points
c  ( tau(i), gtau(i) ), i=1,...,ntau. the b-spline coefficients
c  b c o e f   of the approximating spline are determined from the
c  normal equations using cholesky's method.
c
c******  i n p u t  ******
c  t(1), ..., t(n+k)  the knot sequence
c  n.....the dimension of the space of splines of order k with knots t.
c  k.....the order
c
c  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
c        sed by the arbitrary dimension statement for  biatx  below, but
c        is  n o w h e r e   c h e c k e d   for.
c
c******  w o r k  a r r a y s  ******
c  q....a work array of size (at least) k*n. its first  k  rows are used
c       for the  k  lower diagonals of the gramian matrix  c .
c  diag.....a work array of length  n  used in bchfac .
c
c******  i n p u t  via  c o m m o n  /data/  ******
c  ntau.....number of data points
c  (tau(i),gtau(i)), i=1,...,ntau     are the  ntau  data points to be
c        fitted .
c
c******  o u t p u t  ******
c  bcoef(1), ..., bcoef(n)  the b-spline coeffs. of the penalized l2-appr.
c
c******  m e t h o d  ******
c  the b-spline coefficients of the l2-appr. are determined as the sol-
c  ution of the normal equations
c     sum ( (b(i),b(j))*bcoef(j) : j=1,...,n)  = (b(i),g),
c                                               i = 1, ..., n .
c  here,  b(i)  denotes the i-th b-spline,  g  denotes the function to
c  be approximated, and the  i n n e r   p r o d u c t  of two funct-
c  ions  f  and  g  is given by
c      (f,g)  :=  sum ( f(tau(i))*g(tau(i)) : i=1,...,ntau) .
c  the arrays  t a u  and  w e i g h t  are given in common block
c   d a t a , as is the array  g t a u  containing the sequence
c  g(tau(i)), i=1,...,ntau.
c  the relevant function values of the b-splines  b(i), i=1,...,n, are
c  supplied by the subprogram  b s p l v b .
c     the coeff.matrix  c , with
c           c(i,j)  :=  (b(i), b(j)) + lambda * omega(i, j), i,j=1,...,n,
c  of the normal equations is symmetric and (2*k+1)-banded, therefore
c  can be specified by giving its k+1 bands at or below the diagonal. for
c  i=1,...,n,  we store
c   (b(i),b(j)) + lambda * omega(i, j)  =  c(i,j)  in  q(i-j+1,j), j=max(1, i-k), ... ,i
c  and the right side
c   (b(i), g )  in  bcoef(i) .
c  since b-spline values are most efficiently generated by finding sim-
c  ultaneously the value of  e v e r y  nonzero b-spline at one point,
c  the entries of  c  (i.e., of  q ), are generated by computing, for
c  each ll, all the terms involving  tau(ll)  simultaneously and adding
c  them to all relevant entries, then adding lambda * omega(i, j), computed using 
c  the get_omega function.
*/
      
	  const int kmax = 20;
	  int k, n, i, j, jj, left, leftmk, ll, mm, ntau;
	  double * bcoef, *diag, ** q, *t, biatx[21], dw;
	  double *tau, *gtau;
		        
        n = sp->n;
        k = sp->order;
	  t = sp->knot -1;
	  ntau = data->nobs;
	  tau = data->xdata -1;
	  gtau = data->ydata -1;
	  q = get_mat(k+2, n+1);
	  diag =  (double *) CALLOC(n+1, sizeof(double));
        bcoef = (double *) CALLOC(n+1, sizeof(double));
        
        left = k;
        leftmk = 0;
	  for(ll=1; ll<= ntau; ll++){
      
   /*   locate  l e f t  s.t. tau(ll) in (t(left),t(left+1)) */
label_10:   if (left == n)            goto label_15;
            if (tau[ll] < t[left+1])  goto label_15;
            left = left+1;
            leftmk = leftmk + 1;
			goto label_10;
label_15:   bsplvb ( t, k, 1, tau[ll], left, biatx );
/*
c        biatx(mm) contains the value of b(left-k+mm) at tau(ll).
c        hence, with  dw := biatx(mm), the number dw*gtau(ll)
c        is a summand in the inner product
c           (b(left-k+mm), g)  which goes into  bcoef(left-k+mm)
c        and the number biatx(jj)*dw is a summand in the inner product
c           (b(left-k+jj), b(left-k+mm)), into  q(jj-mm+1,left-k+mm)
c        since  (left-k+jj) - (left-k+mm) + 1  =  jj - mm + 1 .
*/
         
		for(mm=1; mm<=k; mm++){
            dw = biatx[mm];
            j = leftmk + mm;
            bcoef[j] = dw*gtau[ll] + bcoef[j];
            i = 1;
         
				for(jj=mm; jj<=k; jj++){
                  q[i][j] = biatx[jj]*dw + q[i][j];
                  i = i + 1;
				}
			 }
	  }
/*       adjust q to account for the penalty */
/*       note: penalty only applies to piecewise terms, not x^p terms */

        for (i = 1; i<= k + 1; i++)
             for (j = 1; j <= n + 1 - i; j++)
                  q[i][j] += lambda*get_omega(t, n, k, i+j-1, j);

/*
c
c             construct cholesky factorization for  c  in  q , then use
c             it to solve the normal equations
c                    c*x  =  bcoef
c             for  x , and store  x  in  bcoef .
*/
      bchfac ( q, k + 1, n, diag );
      bchslv ( q, k + 1, n, bcoef );
	  for(i=0; i<sp->n; i++) 
		  sp->coef[i] = bcoef[i+1];
	  free_mat(q), free(diag), free(bcoef);
return;
}

void bchinvb( double **w , int nbands, int nrow, double *diag )
{
/* This function efficiently computes the first nrow bands of 
   the inverse of a (symmetric positive definite) banded matrix 
   after a Cholesky decomposition has been performed.  It uses an algorithm
   found in Hutchinson and deHoog (1985). */

int i, j, ell;
double temp;

/* case i = nrow is a special case  - no code needed */

/* temp = w[1][nrow];
w[1][nrow] = temp; */

for (i = nrow - 1; i >= 1; i--)
{
    /* copy elements of L into diag */
    for (ell = 1; ell < nbands; ell++)
        if (ell <= nrow - i)
             diag[ell] = w[ell + 1][i];

    for (j = 1; j < nbands; j++)
        if (j <= nrow - i)
        {

        /* compute element (i, i + j) of the inverse and store it*/
        temp = 0.0;

        for (ell = 1; ell < nbands; ell++)
            if (ell <= nrow - i) 
                {
                if (ell <= j)
                    temp -= diag[ell] * w[j + 1 - ell][i + ell];
                else
                    temp -= diag[ell] * w[ell + 1 - j][i + j];
                }
        w[j + 1][i] = temp;
        }

  
    /* compute element (i, i) of the inverse and store it*/
    temp = w[1][i];

    for (ell = 1; ell < nbands; ell++)
        if (ell <= nrow - i)
             temp -= diag[ell] * w[ell + 1][i];
 
    w[1][i] = temp;
}

}

double trace_hat_matrix( struct L2_1D_DATA * data, struct SP_1D * sp , double lambda)
/* this function computes the trace of the hat matrix
X(X^T X + lambda omega)^{-1} X^T, which is also the trace of 
(X^T X + lambda omega)^{-1} X^T X */
{
  
	  const int kmax = 20;
	  int k, n, i, j, jj, left, leftmk, ll, mm, ntau;
	  double *diag, ** q, ** xtx, *t, biatx[21], dw;
	  double *tau, *gtau;
        double trace = 0.0;

       /* trace = sp->n if no penalty */
       if (lambda == 0.0)
           return sp->n;
		        
        n = sp->n;
        k = sp->order;
	  t = sp->knot -1;
	  ntau = data->nobs;
	  tau = data->xdata -1;
	  gtau = data->ydata -1;
	  q = get_mat(k+2, n+1);
	  xtx = get_mat(k+2, n+1);

	  diag =  (double *) CALLOC(n+1, sizeof(double));
      
        left = k;
        leftmk = 0;
	  for(ll=1; ll<= ntau; ll++){
      
   /*   locate  l e f t  s.t. tau(ll) in (t(left),t(left+1)) */
label_10:   if (left == n)            goto label_15;
            if (tau[ll] < t[left+1])  goto label_15;
            left = left+1;
            leftmk = leftmk + 1;
			goto label_10;
label_15:   bsplvb ( t, k, 1, tau[ll], left, biatx );
/*
c        biatx(mm) contains the value of b(left-k+mm) at tau(ll).
c        hence, with  dw := biatx(mm), the number dw*gtau(ll)
c        is a summand in the inner product
c           (b(left-k+mm), g)  which goes into  bcoef(left-k+mm)
c        and the number biatx(jj)*dw is a summand in the inner product
c           (b(left-k+jj), b(left-k+mm)), into  q(jj-mm+1,left-k+mm)
c        since  (left-k+jj) - (left-k+mm) + 1  =  jj - mm + 1 .
*/
         
		for(mm=1; mm<=k; mm++){
            dw = biatx[mm];
            j = leftmk + mm;
            i = 1;
         
				for(jj=mm; jj<=k; jj++){
                  q[i][j] = biatx[jj]*dw + q[i][j];
                  i = i + 1;
				}
			 }
	  }

       /* store the values of q in xtx so that they are available later */
       for (i = 1; i<= k; i++)
             for (j = 1; j <= n + 1 - i; j++)
                  xtx[i][j] = q[i][j];



/*       adjust q to account for the penalty */
/*       note: penalty only applies to piecewise terms, not x^p terms */

        for (i = 1; i<= k + 1; i++)
             for (j = 1; j <= n + 1 - i; j++)
                  q[i][j] += lambda*get_omega(t, n, k, i+j-1, j);

      /*  construct cholesky factorization for  q  */
      bchfac ( q, k + 1, n, diag );

      /* compute the first k + 1 bands of the inverse of q */
      bchinvb ( q, k + 1, n, diag );

      /* compute trace from q^{-1} and xtx */
      /* i = 1 */
      for (j = 1; j <= n; j++)
           trace += q[1][j] * xtx[1][j];

      for (i = 2; i<= k; i++)
             for (j = 1; j <= n + 1 - i; j++)
                   trace += 2 * q[i][j] * xtx[i][j];
                 
      free_mat(q), free_mat(xtx), free(diag);
      return trace;
}

double trace_hat_matrix_fit( struct L2_1D_DATA * data, struct SP_1D * sp , double lambda)
/* This function computes the trace of the hat matrix
X(X^T X + lambda omega)^{-1} X^T, which is also the trace of 
(X^T X + lambda omega)^{-1} X^T X.  It also fits the spline at the same time, thus avoiding 
computing the penalty matrix twice. */
{
	  const int kmax = 20;
	  int k, n, i, j, jj, left, leftmk, ll, mm, ntau;
	  double *bcoef, *diag, ** q, ** xtx, *t, biatx[21], dw;
	  double *tau, *gtau;
        double trace = 0.0;
		        
        n = sp->n;
        k = sp->order;
	  t = sp->knot -1;
	  ntau = data->nobs;
	  tau = data->xdata -1;
	  gtau = data->ydata -1;
	  q = get_mat(k+2, n+1);
	  xtx = get_mat(k+2, n+1);

	  diag =  (double *) CALLOC(n+1, sizeof(double));
        bcoef = (double *) CALLOC(n+1, sizeof(double));
      
        left = k;
        leftmk = 0;
	  for(ll=1; ll<= ntau; ll++){
      
   /*   locate  l e f t  s.t. tau(ll) in (t(left),t(left+1)) */
label_10:   if (left == n)            goto label_15;
            if (tau[ll] < t[left+1])  goto label_15;
            left = left+1;
            leftmk = leftmk + 1;
			goto label_10;
label_15:   bsplvb ( t, k, 1, tau[ll], left, biatx );
/*
c        biatx(mm) contains the value of b(left-k+mm) at tau(ll).
c        hence, with  dw := biatx(mm), the number dw*gtau(ll)
c        is a summand in the inner product
c           (b(left-k+mm), g)  which goes into  bcoef(left-k+mm)
c        and the number biatx(jj)*dw is a summand in the inner product
c           (b(left-k+jj), b(left-k+mm)), into  q(jj-mm+1,left-k+mm)
c        since  (left-k+jj) - (left-k+mm) + 1  =  jj - mm + 1 .
*/
         
		for(mm=1; mm<=k; mm++){
            dw = biatx[mm];
            j = leftmk + mm;
            bcoef[j] = dw*gtau[ll] + bcoef[j];
            i = 1;
         
				for(jj=mm; jj<=k; jj++){
                  q[i][j] = biatx[jj]*dw + q[i][j];
                  i = i + 1;
				}
			 }
	  }

       /* store the values of q in xtx so that they are available later */ 
        for (i = 1; i<= k; i++)
             for (j = 1; j <= n + 1 - i; j++)
                  xtx[i][j] = q[i][j];



/*       adjust q to account for the penalty */
/*       note: penalty only applies to piecewise terms, not x^p terms */

        for (i = 1; i<= k + 1; i++)
             for (j = 1; j <= n + 1 - i; j++)
                  q[i][j] += lambda*get_omega(t, n, k, i+j-1, j);


/*
c
c             construct cholesky factorization for  c  in  q , then use
c             it to solve the normal equations
c                    c*x  =  bcoef
c             for  x , and store  x  in  bcoef .
*/
      bchfac ( q, k + 1, n, diag );
      bchslv ( q, k + 1, n, bcoef );
	  for(i=0; i<sp->n; i++) 
		  sp->coef[i] = bcoef[i+1];


      /* compute the first k + 1 bands of the inverse of q */
      bchinvb ( q, k + 1, n, diag );

      /* compute trace from q^{-1} and xtx */
      /* i = 1 */
      for (j = 1; j <= n; j++)
           trace += q[1][j] * xtx[1][j];

      for (i = 2; i<= k; i++)
             for (j = 1; j <= n + 1 - i; j++)
                   trace += 2 * q[i][j] * xtx[i][j];
                 
      free_mat(q), free_mat(xtx), free(diag), free(bcoef);
      return trace;
}


double get_L2_error(struct SP_1D * sp, struct L2_1D_DATA * data ){
int i;
double sse;	
sse = 0.;
for(i=0; i<data->nobs; i++)
	sse += pow(data->ydata[i] - sp_1d_value(sp, data->xdata[i], 0), 2);

return sse;
}

double gcv( struct L2_1D_DATA * data, struct SP_1D * sp, double lambda, int minspace)
{
    double sse;
    double trace;
    double n;
    double answer;

    /* Set GCV = INFINITY if knots too close together */
    if (check_knots(data, sp, minspace))
         return INFINITY;

    sse = get_L2_error(sp, data);
    trace = trace_hat_matrix(data, sp, lambda);
    n = (double) data->nobs;
    answer = (sse/n) / (pow( (n - trace)/n, 2.0));

    return answer;
}

double gcv_fit( struct L2_1D_DATA * data, struct SP_1D * sp, double lambda, int minspace)
{
/* computes GCV and fits spline at the same time */
    double sse;
    double trace;
    double n;
    double answer;

    /* Set GCV = INFINITY if knots too close together */
    if (check_knots(data, sp, minspace))
         return INFINITY;

    trace = trace_hat_matrix_fit(data, sp, lambda);
    sse = get_L2_error(sp, data);
    n = (double) data->nobs;
    answer = (sse/n) / (pow( (n - trace)/n, 2.0));

    return answer;
}


/* Checks spacing between knots
If spacing too small, return 1
Otherwise, return 0 */

int check_knots(struct L2_1D_DATA * data, struct SP_1D * sp, int minspace)
{
    int i, j, curknot, count;

    j = 0;

    for (i = 0; i <= sp->n - sp->order; i++)
    {
           count = 0;
           curknot = sp->order + i - 1;
           while (data->xdata[j] < sp->knot[curknot + 1])
           {
                j++;
                count++;
           }
           if (count < minspace)
                return 1;
    }
    return 0;
}

void swap_sp(struct SP_1D ** sp1, struct SP_1D ** sp2)
{
     struct SP_1D * temp;

     temp = *sp1;
     *sp1 = *sp2;
     *sp2 = temp;
}

void swap_scalar(double* x1, double* x2)
{
     double temp;

     temp = *x1;
     *x1 = *x2;
     *x2 = temp;
}
     
int qsort_data_partition(struct L2_1D_DATA * out, int left, int right, int pivotIndex)
{
     int storeIndex, i;
     double pivotValue;
     
     pivotValue = out->xdata[pivotIndex];

     /* Move pivot to end */
     swap_scalar(&(out->xdata[pivotIndex]), &(out->xdata[right]));
     swap_scalar(&(out->ydata[pivotIndex]), &(out->ydata[right]));

     storeIndex = left - 1;

     for(i = left; i <= right-1; i++)
         if (out->xdata[i] <= pivotValue)
         {
             storeIndex++;
             swap_scalar(&(out->xdata[storeIndex]), &(out->xdata[i]));
             swap_scalar(&(out->ydata[storeIndex]), &(out->ydata[i]));
         }

     /* Move pivot to its final place */
     swap_scalar(&(out->xdata[right]), &(out->xdata[storeIndex + 1]));
     swap_scalar(&(out->ydata[right]), &(out->ydata[storeIndex + 1]));

     return (storeIndex + 1);
}

void qsort_data(struct L2_1D_DATA * out, int left, int right)
{
     int pivotIndex, pivotNewIndex;

     if (right > left)
     {
           pivotIndex = left;
           pivotNewIndex = qsort_data_partition(out, left, right, pivotIndex);
           qsort_data(out, left, pivotNewIndex - 1);
           qsort_data(out, pivotNewIndex + 1, right);
     }
}


void sort_data(struct L2_1D_DATA * in, struct L2_1D_DATA * out)
{
    int i, nobs, sort_flag;

    nobs = in->nobs;

    for(i = 0; i < nobs; i++)
    {
        out->xdata[i] = in->xdata[i];
        out->ydata[i] = in->ydata[i];
    }    

    /* Check to see if data sorted */
    sort_flag = 1; 

    for (i = 0; i < nobs - 1; i++)
        if (out->xdata[i] > out->xdata[i + 1])
        { 
           sort_flag = 0;
           break;
        }

    if (sort_flag == 0)
        qsort_data(out, 0, nobs - 1);
}

void scale_data(struct L2_1D_DATA * in, struct L2_1D_DATA * out, double * shift, double * scale)
{
     int i, nobs;
     double lowx, highx;

     nobs = in->nobs;
 
     lowx = in->xdata[0];
     highx = in->xdata[0];

     for(i = 1; i<nobs; i++)
     {
         if (in->xdata[i] < lowx) lowx = in->xdata[i];
         if (in->xdata[i] > highx) highx = in->xdata[i];
     }

     *shift = lowx;

     *scale = highx - lowx;

     for(i = 0; i<nobs; i++)
     {
         out->xdata[i] = (in->xdata[i] - *shift)/(*scale);
         out->ydata[i] = in->ydata[i];
     }    
}
     
double GSJS(struct L2_1D_DATA * data)
{
    int i, n = data->nobs;
    double A, B, epsilon, answer = 0.0;

    for (i = 1; i < n-1; i++)
    {
          A = (data->xdata[i + 1] - data->xdata[i])/
              (data->xdata[i + 1] - data->xdata[i - 1]);
          B = 1.0 - A;
          epsilon = (data->ydata[i] - A*data->ydata[i-1] - 
                    B*data->ydata[i+1])/(pow(1 + A*A + B*B, 0.5));
          answer += epsilon * epsilon;
    }
    
    answer /= ((double) n - 2.0);

    return answer;
}

double get_chg_basis_mat_order_2(double* t, int i, int j)
{
      /* This function only works for the order 2 (linear) case, and if
     j >= 2 (bottom portion of matrix) */

      if (i <= j - 3)
          return 0.0;
      if (i == j - 2)
          return 1/(t[i+2] - t[i+1]);
      if (i == j - 1)
          return - 1/(t[i+1] - t[i]) - 1/(t[i+2] - t[i+1]);
      if (i == j)
          return 1/(t[i+1] - t[i]);
      if (i > j)
          return 0.0;
      return 0.0;      
}


double get_chg_basis_mat_order_3(double* t, int i, int j)
{
      /* This function only works for the order 3 (quadratic) case, and
      if j >= 3 (bottom portion of matrix) */

      if (i <= j - 4)
          return 0.0;
      if (i == j - 3)
          return -1/((t[i+3] - t[i+1]) * (t[i+3] - t[i+2]));
      if (i == j - 2)
          return (1/(t[i+3]-t[i+2]) + 1/(t[i+2]-t[i]))/(t[i+2]-t[i+1]);
      if (i == j - 1)
          return (1/(t[i+1]-t[i+3]) - 1/(t[i+1]-t[i]))/(t[i+2]-t[i+1]);
      if (i == j)
          return  1/((t[i]-t[i+1]) * (t[i]-t[i+2]));
      if (i > j)
          return 0.0;
      return 0.0;      
}


double get_chg_basis_mat_order_4(double* t, int i, int j)
{
      /* This function only works for the order 4 (cubic) case, and if
       j >= 4  (bottom portion of matrix) */

      if (i <= j - 5)
          return 0.0;
      if (i == j - 4)
          return 1/((t[i+4]-t[i+1]) * (t[i+4]-t[i+2]) * 
                 (t[i+4]-t[i+3]));
      if (i == j - 3)
          return (-1/(t[i+4]-t[i+3]) + 
                 1/(t[i]-t[i+3]))/((t[i+3]-t[i+2]) *(t[i+3]-t[i+1]));
      if (i == j - 2)
          return (1/(t[i+4]-t[i+2]) - 
                  1/(t[i]-t[i+2]))/((t[i+2]-t[i+1]) * (t[i+3]-t[i+2]));
      if (i == j - 1)
          return (-1/(t[i+4]-t[i+1]) + 
                 1/(t[i]-t[i+1]))/((t[i+2]-t[i+1]) * (t[i+3]-t[i+1]));
      if (i == j)
          return  -1/((t[i]-t[i+1]) * (t[i]-t[i+2]) * (t[i]-t[i+3]));
      if (i > j)
          return 0.0;      
      return 0.0;
}
