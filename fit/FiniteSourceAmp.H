#include <math.h>
#include <stdlib.h>
#define PI 3.1415927
#define EPSM 1.0e-3
#define JMAX 13


/******************************************************

FinSrcAmp.h: set of routines for calculating microlensing
		 amplitude for finite size source stars. Only function
		 needed by user is FS_Amp.

*******************************************************/

double mag(double u, double r);	/* mag(u) for const. surface brightness star */

double mag_noapprox(double u, double r);  /* mag from W&M eq. 9 or 11 */

double mag_large_u(double u, double r); /* mag from sum in eq A4 in Witt&Mao */

double LDmag_large_u(double u, double r, double limba, double limbb); /* limb drk. approx*/

double mag_small_u(double u, double r); /* mag from sum in eq A7 in Witt&Mao */

double mag_func(double u, double r, double u_star, double limba, double limbb);
	/* derivitive of mag. (dm/dr) of const surface brightness star,
	   integrated by parts -- this is integrated over star surface
	   to give mag with limb darkening */

void donothing(void);   /* used to fool optimizing compiler, seperates
			   two commands that can't be combined */

double midpntz(double u, double u_star, double limba, double limbb, int n);
	/* modified N.R. integration routine */

double FS_Amp(double u, double u_star, double limba, double limbb);
	/* modified N.R. integration routine -- returns amplitude
	   from finite source star microlensing */

/*************************************************************************

NOTE ---  ellpi(PI/2, en, ak) = ellpip(en, ak)   =   ellpisum(-en, ak)
	  [Numerical Recipes]   [Numerical Recipes]  [Gradshteyn & Ryzhik]

	  because of different definitions of ellpi


	 N.R.:  ellpi = int(dx [1 + n*sin(x)^2]^(-1) [1  - (k sin(x))^2]^(-1/2)
	 G.R.:  ellpi = int(dx [1 - n*sin(x)^2]^(-1) [1  - (k sin(x))^2]^(-1/2)

**************************************************************************/
#ifndef ELLIPTIC
#define ELLIPTIC
						/* elliptic integral of:  */
double ellf(double phi, double ak);			/* first kind     */
double elle(double phi, double ak);			/* second kind    */
double ellpi(double phi, double en, double ak);		/* third kind     */


					/* perfect elliptic integral of:  */
double ellfp(double ak);				/* first kind     */
double ellep(double ak);				/* second kind    */
double ellpip(double en, double ak);			/* third kind     */


			/* summation formula for perfect integral of:     */
double ellfsum(double ak);				/* first kind     */
double ellesum(double ak);				/* second kind    */
double ellpisum(double en, double ak);			/* third kind     */


		      /* factorial-type functions for summation formulas  */
unsigned long oddfacd(int m);				/* (2m - 1)!!     */
unsigned long evenfacd(int m);				/* (2m)!!         */
unsigned long fac(int m);				/* m!             */

#endif

#include <stdio.h>
#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0
 
double FS_dAmp (double (*func)(double, double, double, double), int which,
		double var0, double var1, double var2, double var3,
		double h, double *err);
        /* Modified dfridr derivative routine */

double **dmatrix(int nrl,int nrh,int ncl,int nch);

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
