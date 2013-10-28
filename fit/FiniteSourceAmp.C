#include <math.h>
#include <stdio.h>
#include "FiniteSourceAmp.H"

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))

static double minarg1, minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
(minarg1) : (minarg2))

double mag(double u, double r)
{
   if(r > 32) {
      if(u > r) return mag_large_u(u, r);
      else return mag_small_u(u, r);
   } else {
      if(u > 0.90625*r + 3.0) return mag_large_u(u, r);
      else if(r > 0.90625*u + 3.0) return mag_small_u(u, r);
      else return mag_noapprox(u, r);
   }
}


double mag_noapprox(double u, double r)
{
   double c1, c2, c3, fac, k, n, mu;

   n = 4*u*r / ((u+r) * (u+r));
   if(n < 0.999999) {

      fac = sqrt(4 + (u-r)*(u-r));

      c1 = (u+r) / (2*r*r) * fac;
      c2 = (u-r) / (r*r) * (4 + (u*u - r*r)/2) / fac;
      c3 = 2*(u-r)*(u-r) / (r*r*(u+r)) * (1 + r*r) / fac;

      k = 2*sqrt(n)/fac;
      if (n < 0)
	 printf ("Magnoappx, n = %lf, r = %lf, u = %lf\n", n, r, u);

      mu = (c1*ellep(k)-c2*ellfp(k)+c3*ellpip(-n, k))/PI;

   } else {

      mu = (2/r + (1+r*r)/(r*r)*(PI/2 + asin((r*r - 1)/(r*r +1))))/PI;
   
   }

   return mu;
}




double LDmag_large_u(double u, double r, double limba, double limbb)
{
   double c0, c2, c4, mu, denom;
   double u2, u3, u4, u5, u6;

   u2 = u*u;
   u3 = u*u2;
   u4 = u*u3;
   u5 = u*u4;
   u6 = u*u5;

   c0 = (u2 + 2.0) / u / sqrt(u2 + 4.0);
   c2 = 4.0 * (u2 + 1.0) / u3 / pow(u2 + 4.0, 2.5);
   c4 = 2.0 * (12.0 + 14.0*u2 + 6.0*u4 + 3.0*u6) / u5 / pow(u2 + 4.0, 4.5);

   if (limbb == 0.0){
      c2 *= (15.0 - 7.0*limba) / (15.0 - 5.0*limba);
      c4 *= (105.0 - 57.0*limba) / (105.0 - 35.0*limba);
   }
   else {
      denom = -6 + 2*limba + limbb;
      c2 *= (-6 + 14/5*limba + 8/5*limbb) / denom;
      c4 *= (-6 + 114/35*limba + 141/70*limbb) / denom;
   }
   
   mu = c0 + r*r*(c2 + c4*r*r);

   return mu;
}


double mag_large_u(double u, double r)
{
   double c0, c2, c4, mu;
   double u2, u3, u4, u5, u6;

   u2 = u*u;
   u3 = u*u2;
   u4 = u*u3;
   u5 = u*u4;
   u6 = u*u5;

   c0 = (u2 + 2.0) / u / sqrt(u2 + 4.0);
   c2 = 4.0 * (u2 + 1.0) / u3 / pow(u2 + 4.0, 2.5);
   c4 = 2.0 * (12.0 + 14.0*u2 + 6.0*u4 + 3.0*u6) / u5 / pow(u2 + 4.0, 4.5);

   mu = c0 + r*r*(c2 + c4*r*r);

   return mu;
}

double mag_small_u(double u, double r)
{
   double a0, a2, a4, mu;
   double r2, r3, r4, r5, fac;

   r2 = r*r;
   r3 = r2*r;
   r4 = r3*r;
   r5 = r4*r;

   fac = sqrt(r2 + 4.0);

   a0 = fac / r;
   a2 = -4.0 / r3 / pow(fac, 3.0);
   a4 = -6.0 * (r4 + 2.0*r2 + 2.0) / r5 / pow(fac, 7.0);

   mu = a0 + u*u*(a2 + a4*u*u);

   return mu;
}

double mag_func(double u, double r, double u_star, double limba, double limbb)
{
   double h, dudr, temp;
   void donothing(void);
   double mag(double u, double r);

   h = FMAX(0.005*r, 1.0e-7);

   /*** force h to be representable number ***/

   temp = r + h;
   donothing();
   h = temp - r;

   if(r-h >0) dudr = (mag(u, r+h) - mag(u, r-h)) / 2.0 / h;
   else dudr = (mag(u, r+h) - mag(u, r)) / h;
   /*
   printf("%f %f %f %f %f\n", u, r, u_star, h, dudr);
   */

/*      return r*(mag(u,r) + r/2.0 * dudr) *
	(1.0 - limbd + limbd*sqrt(1.0 - r*r/u_star/u_star)); */
   
   return r*(mag(u,r) + r/2.0 * dudr) *
	   (1.0 - limba*(1.0 - sqrt(1.0 - r*r/u_star/u_star)) -
	    limbb*(1.0 - sqrt(1.0 - r*r/u_star/u_star))*
	    (1.0 - sqrt(1.0 - r*r/u_star/u_star)));
}


void donothing(void){return;}


double midpntz(double u, double u_star, double limba, double limbb, int n)
{
   double r, tnm, sum, del, ddel;
   static double s;
   int it, j;

   if(n == 1) return (s =u_star*mag_func(u, u_star/2.0, u_star, limba, limbb));
   else {
      for(it=1, j=1; j < n-1; j++) it *= 3;
      tnm = it;
      del = u_star/(3.0*tnm);
      ddel = del + del;
      r = 0.5*del;
      sum = 0.0;
      for(j = 1; j <= it; j++)
      {
	 sum += mag_func(u, r, u_star, limba, limbb);
	 r += ddel;
	 sum += mag_func(u, r, u_star, limba, limbb);
	 r += del;
      }
      s = (s + u_star*sum/tnm) / 3.0;
      return s;
   }
}

double FS_Amp(double u, double u_star, double limba, double limbb)
{
   double midpntz(double u, double u_star, double limba, double limbb, int n);
   int j;
   double s, st, ost, os;
   double r_eff;
   double D;

   if (u < 0)
      printf ("u = %f\n", u);

   /* r_eff = u_star * sqrt(1.0 - limbd/3.0); */
   r_eff = u_star * sqrt(1.0 - limba/3.0 - limbb/6.0);
   
   if(r_eff > 50.0) return 1.0;   /* max amp < 1.001 */

   if( (limba == 0) && (limbb == 0) ) return mag(u, u_star);
                                  /* uniform surface function */

   if(u > u_star+3.0) return LDmag_large_u(u, u_star, limba, limbb);

   /* D = u_star*u_star*(.5 - limbd/6.0); */
   D = u_star*u_star*(.5 - limba/6.0 - limbb/12.0);
      
   ost = os = -1.0e30;
   for(j = 1; j <= JMAX; j++)
   {
      st = midpntz(u, u_star, limba, limbb, j);
      s = (9.0*st - ost)/8.0;
      if(fabs(s - os) < EPSM*fabs(os)) return s/D;
      if(s == 0.0 && os == 0.0 && j > 6) return s;
      os = s;
      ost = st;
   }

   fprintf (stderr, "Too many steps in routine Amp");
   return 1;

/* nrerror("Too many steps in routine Amp");
   return 0.0; */
}

/******************************************************

old version with Simpson's rule


double FS_Amp(double u, double u_star, double limbd)
{
   double midpntz(double u, double u_star, double limbd, int n);
   int j;
   double s, olds;
   double D;

   if(limbd == 0) return mag(u, u_star);

   D = u_star*u_star*(.5 - limbd/6.0);
   olds = -1.0e30;

   for(j = 1; j <= JMAX; j++) {
      s = midpntz(u, u_star, limbd, j);
      if(fabs(s - olds) < EPSM*fabs(olds)) return s;
      if(s == 0.0 && olds == 0.0 && j > 6) return s;
      olds = s;
   }

   nrerror("Too many steps in routine Amp");
   return 0.0;
}

*****************************************************************/
/*******************************************************************

elliptic.c --- library of functions to calculate elliptic
	       integrals of the first, second and third kind
	       to 7 significant digits

	       from Numerical Recipes in C - Second edition

   1st kind: F(phi, k) = double ellf(double phi, double ak)
   2nd kind: E(phi, k) = double elle(double phi, double ak)
   3rd kind: PI(phi, n, k) = double ellpi(double phi, double en, double ak)


			 Matt Lehner     3/29/95

       modified to do perfect integrals     MJL 3/29/95

       **NOTE**  ellpi(phi, en, ak) in N.R. = PI(phi, -n, k) in G.R.
		  MJL 4/11/95

********************************************************************/

#define PI 3.1415927

/************************ constants for rc *************************/

#define ERRTOLC 0.04
#define TINYC 1.69e-38
#define SQRTNY 1.3e-19
#define BIGC 3.0e37
#define TNBG (TINYC*BIGC)
#define COMP1 (2.236/SQRTNY)
#define COMP2 (TNBG*TNBG/25.0)
#define THIRD (1.0/3.0)
#define C1C 0.3
#define C2C (1.0/7.0)
#define C3C 0.375
#define C4C (9.0/22.0)


/************************ constants for rd *************************/

#define ERRTOLD 0.05
#define TINYD 1.0e-25
#define BIGD 4.5e21
#define C1D (3.0/14.0)
#define C2D (1.0/6.0)
#define C3D (9.0/22.0)
#define C4D (3.0/26.0)
#define C5D (0.25*C3D)
#define C6D (1.5*C4D)


/************************ constants for rf *************************/

#define ERRTOLF 0.08
#define TINYF 1.5e-38
#define BIGF 3.0e37
#define C1F (1.0/24.0)
#define C2F 0.1
#define C3F (3.0/44.0)
#define C4F (1.0/14.0)


/************************ constants for rj *************************/

#define ERRTOLJ 0.05
#define TINYJ 2.5e-13
#define BIGJ 9.0e11
#define C1J (3.0/14.0)
#define C2J (1.0/3.0)
#define C3J (3.0/22.0)
#define C4J (3.0/26.0)
#define C5J (0.75*C3J)
#define C6J (1.5*C4J)
#define C7J (0.5*C2J)
#define C8J (C3J+C3J)

double rc(double x, double y)
{
   double alamb, ave, s, w, xt, yt;

   if (x < 0.0 || y == 0.0 || (x+fabs(y)) < TINYC || (x+fabs(y)) > BIGC ||
       (y < -COMP1 && x > 0.0 && x < COMP2)) {
	   fprintf(stderr, "%.10f %.10f\n", x, y);
	   /*nrerror("invalid arguments in rc"); */
           return 1;
   }

   if(y > 0.0)
   {
      xt = x;
      yt = y;
      w = 1.0;
   }
   else
   {
      xt = x-y;
      yt = -y;
      w = sqrt(x)/sqrt(xt);
   }

   do {
      alamb = 2.0*sqrt(xt)*sqrt(yt) + yt;
      xt = 0.25*(xt+alamb);
      yt = 0.25*(yt+alamb);
      ave = THIRD*(xt+yt+yt);
      s = (yt - ave)/ave;
   } while (fabs(s) > ERRTOLC);

   return w*(1.0+s*s*(C1C+s*(C2C+s*(C3C+s*C4C))))/sqrt(ave);
}



double rd(double x, double y, double z)
{
   double alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty,
          sqrtz, sum, xt, yt, zt;

   if (FMIN(x,y) < 0 || FMIN(x+y, z) < TINYD || FMAX(FMAX(x,y),z) > BIGD)
	   /*nrerror("invalid arguments in rd");*/
           return 1;
   
   xt = x;
   yt = y;
   zt = z;
   sum = 0.0;
   fac = 1.0;

   do {
      sqrtx = sqrt(xt);
      sqrty = sqrt(yt);
      sqrtz = sqrt(zt);
      alamb = sqrtx*(sqrty+sqrtz) + sqrty*sqrtz;
      sum += fac/(sqrtz * (zt+alamb));
      fac = 0.25*fac;
      xt = 0.25*(xt+alamb);
      yt = 0.25*(yt+alamb);
      zt = 0.25*(zt+alamb);
      ave = 0.2*(xt+yt+3.0*zt);
      delx = (ave - xt)/ave;
      dely = (ave - yt)/ave;
      delz = (ave - zt)/ave;
   } while (FMAX(FMAX(fabs(delx),fabs(dely)), fabs(delz)) > ERRTOLD);

   ea = delx*dely;
   eb = delz*delz;
   ec = ea - eb;
   ed = ea - 6.0*eb;
   ee = ed + ec + ec;

   return 3.0*sum+fac*(1.0+ed*(-C1D+C5D*ed-C6D*delz*ee)
	  +delz*(C2D*ee+delz*(-C3D*ec+delz*C4D*ea)))/(ave*sqrt(ave));
}
 


double rf(double x, double y, double z)
{
   double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt;

   if (FMIN(FMIN(x,y),z) < 0 || FMIN(FMIN(x+y, x+z), y+z) < TINYF ||
       FMAX(FMAX(x,y),z) > BIGF)
	   /*nrerror("invalid arguments in rf");*/
           return 1;
   
   xt = x;
   yt = y;
   zt = z;

   do {
      sqrtx = sqrt(xt);
      sqrty = sqrt(yt);
      sqrtz = sqrt(zt);
      alamb = sqrtx*(sqrty+sqrtz) + sqrty*sqrtz;
      xt = 0.25*(xt+alamb);
      yt = 0.25*(yt+alamb);
      zt = 0.25*(zt+alamb);
      ave = THIRD*(xt+yt+zt);
      delx = (ave - xt)/ave;
      dely = (ave - yt)/ave;
      delz = (ave - zt)/ave;
   } while (FMAX(FMAX(fabs(delx),fabs(dely)), fabs(delz)) > ERRTOLF);

   e2 = delx*dely - delz*delz;
   e3 = delx*dely*delz;
   return (1.0 + (C1F*e2 - C2F - C3F*e3)*e2 + C4F*e3)/sqrt(ave);
}



double rj(double x, double y, double z, double p)
{
   double rc(double x, double y);
   double rf(double x, double y, double z);

   double a, alamb, alpha, ans, ave, b, beta, delp, delx, dely, delz, ea, eb,
	  ec, ed, ee, fac, pt, rcx, rho, sqrtx, sqrty, sqrtz, sum, tau, xt, yt,
	  zt;

   if (FMIN(FMIN(x,y), z) < 0 || FMIN(FMIN(x+y, x+z),FMIN(y+z, fabs(p))) < TINYJ
       || FMAX(FMAX(x,y), FMAX(z, fabs(p))) > BIGJ)
	   /*nrerror("invalid arguments in rj");*/
           return 1;
   
   sum = 0.0;
   fac = 1.0;
   if(p > 0.0)
   {
      xt = x;
      yt = y;
      zt = z;
      pt = p;
   }
   else
   {
      xt = FMIN(FMIN(x,y), z);
      zt = FMAX(FMAX(x,y), z);
      yt = x+y+z-xt-zt;
      a = 1.0/(yt-p);
      b = a*(zt - yt)*(yt - xt);
      pt = yt+b;
      rho = xt*zt/yt;
      tau = p*pt/yt;
      rcx = rc(rho, tau);
   }

   do {
      sqrtx = sqrt(xt);
      sqrty = sqrt(yt);
      sqrtz = sqrt(zt);
      alamb = sqrtx*(sqrty+sqrtz) + sqrty*sqrtz;
      alpha = SQR(pt*(sqrtx+sqrty+sqrtz) + sqrtx*sqrty*sqrtz);
      beta = pt*SQR(pt+alamb);
      sum += fac*rc(alpha, beta);
      fac = 0.25*fac;
      xt = 0.25*(xt+alamb);
      yt = 0.25*(yt+alamb);
      zt = 0.25*(zt+alamb);
      pt = 0.25*(pt+alamb);
      ave = 0.2*(xt+yt+zt+pt+pt);
      delx = (ave - xt)/ave;
      dely = (ave - yt)/ave;
      delz = (ave - zt)/ave;
      delp = (ave - pt)/ave;
   } while (FMAX(FMAX(fabs(delx),fabs(dely)),
	    FMAX(fabs(delz), fabs(delp))) > ERRTOLJ);

   ea = delx*(dely + delz) + dely*delz;
   eb = delx*dely*delz;
   ec = delp*delp;
   ed = ea - 3.0*ec;
   ee = eb + 2.0*delp*(ea - ec);
   ans = 3.0*sum+fac*(1.0+ed*(-C1J+C5J*ed-C6J*ee)+eb*(C7J+delp*(-C8J+delp*C4J))
	 +delp*ea*(C2J-delp*C3J)-C2J*delp*ec)/(ave*sqrt(ave));
   if(p <= 0.0) ans = a*(b*ans+3.0*(rcx-rf(xt, yt, zt)));

   return ans;
}


/************************Legendre elliptic integrals**********************/


double ellf(double phi, double ak)
{
   double rf(double x, double y, double z);
   double s;

   s = sin(phi);
   return s*rf(SQR(cos(phi)), (1.0 - s*ak)*(1.0 + s*ak), 1.0);
}


double elle(double phi, double ak)
{
   double rd(double x, double y, double z);
   double rf(double x, double y, double z);
   double cc, q, s;

   s = sin(phi);
   cc = SQR(cos(phi));
   q = (1.0 - s*ak)*(1.0 + s*ak);

   return s*(rf(cc, q, 1.0) - (SQR(s*ak))*rd(cc, q, 1.0)/3.0);
}


double ellpi(double phi, double en, double ak)
{
   double rf(double x, double y, double z);
   double rj(double x, double y, double z, double p);

   double cc, enss, q, s;

   s = sin(phi);
   enss = en*s*s;
   cc = SQR(cos(phi));
   q = (1.0 - s*ak)*(1.0 + s*ak);
   return s*(rf(cc, q, 1.0) - enss*rj(cc, q, 1.0, 1.0+enss)/3.0);
}


/**********************perfect Legendre integrals ****************/


double ellfp(double ak)
{
   double rf(double x, double y, double z);
   return rf(0, 1.0 - ak*ak, 1.0);
}


double ellep(double ak)
{
   double rd(double x, double y, double z);
   double rf(double x, double y, double z);

   return (rf(0, 1.0 - ak*ak, 1.0) - ak*ak*rd(0, 1.0 - ak*ak, 1.0)/3.0);
}


double ellpip(double en, double ak)
{
   double rf(double x, double y, double z);
   double rj(double x, double y, double z, double p);

   return rf(0, 1.0 - ak*ak, 1.0) - en*rj(0, 1.0 - ak*ak, 1.0, 1.0+en)/3.0;
}


/*********************** summation formulas for perfect integrals **********
************************ ellpisum(n, k) = ellpip(-n, k)  MJL 4/11/95 ******/

double ellfsum(double ak)
{
   int m;
   double x, y, p, sum, ans;

   sum = 0;

   for(m=1; m<10; m++)
   {
      x = (1.0 * oddfacd(m))/evenfacd(m);
      p = 2.0*m;
      y = x*x*pow(ak,p);
      sum += y;
   }

   ans = PI*(1+sum)/2.0;
   return ans;
}


double ellesum(double ak)
{
   int m;
   double x, y, p, sum, ans;

   sum = 0;

   for(m=1; m<10; m++)
   {
      x = (1.0 * oddfacd(m))/evenfacd(m);
      p = 2.0*m;
      y = x*x*pow(ak,p)/(2*m-1);
      sum += y;
   }

   ans = PI*(1-sum)/2.0;
   return ans;
}
   


double ellpisum(double en, double ak)
{
   int m, j;
   double y, p1, p2, p3, p4, sum, ans, f1, f2, f3, f4;

   sum = 0;

   for(m=0; m<12; m++) for(j=0; j<=m; j++)
   {
      f1 = fac(2*m);
      f2 = fac(2*j);
      f3 = fac(m);
      f4 = fac(j);
      p1 = 2*j;
      p2 = m-j;
      p3 = m;
      p4 = j;

      y = f1*f2*pow(ak,p1)*pow(en,p2)/pow(4.0,p3)/pow(4.0,p4)/f3/f3/f4/f4;
      sum += y;
   }

   ans = PI*sum/2.0;
   return ans;
}

/*************** factorial-type functions for summation formulas ********/


unsigned long oddfacd(int m)
{
   if(m==1) return 1;
   else return (2*m - 1)*oddfacd(m-1);
}


unsigned long evenfacd(int m)
{
   if(m==1) return 2;
   else return 2*m*evenfacd(m-1);
}


unsigned long fac(int m)
{
   if(m==0) return 1;
   else return m*fac(m-1);
}

/*****************************************************************************

  Version of dfridr found in Numerical Recipies, to find derivatives of
  functions, called FS_dAmp.  Adapted by Andy Becker for FS_Amp.

  ***************************************************************************/

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
        int i;
        double **m;
 
        m=(double **)malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m) {
	   fprintf (stderr, "allocation failure 1 in dmatrix()");
	   return 0;
	}
        m -= nrl;
 
        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *)malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i]) {
		   fprintf (stderr, ("allocation failure 2 in dmatrix()"));
		   return 0;
		}
                m[i] -= ncl;
        }
        return m;
}


void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
        int i;
 
        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}


double FS_dAmp (double (*func)(double, double, double, double), int which,
		double var0, double var1, double var2, double var3,
		double h, double *err)
{
   int i,j;
   double errt, fac, hh, **a, ans;

   if (h == 0.0) {
      fprintf (stderr, "h must be nonzero in dfridr");
      return 1;
   }
   
   a=dmatrix(1,NTAB,1,NTAB);
   /* Keeps the value of var0, var1, var2 sent to (*func) > 0 */
   if (which == 0)
      hh = h > var0 ? var0 - 0.00001 : h;
   else if (which == 1)
      hh = h > var1 ? var1 - 0.00001 : h;
   else if (which == 2)
      hh = h > var2 ? var2 - 0.00001: h;
   else {
      fprintf (stderr, "Invalid derivative variable in FS_dAmp");
      return 0;
      }

   if (which == 0)
      a[1][1] = ( (*func)(var0+hh, var1, var2, var3)-(*func)(var0-hh, var1, var2, var3) )/(2.0*hh);
   else if (which == 1)
      a[1][1] = ( (*func)(var0, var1+hh, var2, var3)-(*func)(var0, var1-hh, var2, var3) )/(2.0*hh);

/*   else a[1][1] = ( (*func)(var0, var1, var2+hh)-(*func)(var0, var1, var2-hh) )/(2.0*hh); */
      
   *err = BIG;
   for (i=2; i<=NTAB; i++) {
      hh /= CON;
      if (which == 0)
         a[1][i] = ( (*func)(var0+hh, var1, var2, var3)-(*func)(var0-hh, var1, var2, var3))/(2.0*hh);
      else if (which == 1)
	 a[1][i] = ( (*func)(var0, var1+hh, var2, var3)-(*func)(var0, var1-hh, var2, var3))/(2.0*hh);

/*      else a[1][i] = ( (*func)(var0, var1, var2+hh)-(*func)(var0, var1, var2-hh) )/(2.0*hh); */

      fac = CON2;
      for (j=2; j<=i; j++) {
	 a[j][i] = (a[j-1][i]*fac - a[j-1][i-1])/(fac-1.0);
	 fac *= CON2;
	 errt = FMAX(fabs(a[j][i]-a[j-1][i]), fabs(a[j][i]-a[j-1][i-1]));
	 if (errt <= *err) {
	    *err = errt;
	    ans = a[j][i];
	 }
      }
      if (fabs(a[i][i] - a[i-1][i-1]) >= SAFE*(*err)) {
	 free_dmatrix(a, 1, NTAB, 1, NTAB);
	 return ans;
      }
   }
   free_dmatrix(a, 1, NTAB, 1, NTAB);
   return ans;
}
