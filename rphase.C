/* =============================
   Program Description : phase.c  
   =============================
   [James Reimann,  August 1994]

   Converted into a C++ object by John Doug Reynolds, July 1996
      Removed a potentially serious fortran interface error in super().
      Attempted to make safe for multi-threaded applications.
      Repackaged for use with the MACHO Common Lightcurve Interface.


Function:
---------

This program fits periodic curves to light-curve data, by fitting cosines,
supersmoother, or periodic cubic splines.  It first searches for periods
that minimize a measure of fit (MOF) over a group of periods that is
equally spaced in frequency (the initial grid).  The MOF used is generally
weighted residual sum-of-squares, but weighted sum of absolute residuals is
used for supersmoother.  The program takes the periods corresponding to the
best local mimima of the MOF, and places grids around these to refine the
estimates (the final grids).  The program ignores minima that are close to
a better minimum, so that each estimate should be `independent'.
The program may additionally use the full data in a second refining.
The program returns the best frequencies/periods in order of quality of fit. 


Precision:
----------

For the settings: RATE=6 SUBN=100 FNUM=100 FWID=1,
the frequency estimates are correct to about 4 decimal places.
The period estimates have precision that varies with size: 5sf for 0.04
days, 4sf for 0.4 days, 3sf for 4 days and 2sf for 40 days. This is because
we get more information from the data about short periods than we do about
long periods.


Parameters:
-----------

TYPE:  Type of curve:   1 = trigonometric polynomial
			2 = supersmoother
			3 = cubic spline
HARM:  If TYPE=1, The number of harmonics in the trigonometric model
		  Takes values 1,2,3,...                       
BASS:  If TYPE=2, Bass parameter in supersmoother - range is [0,10]   
                  0 recommended.  
KNOT:  If TYPE=3, Number of knots in cubic spline - need at least 8.
SUBN:  Number of observations in subset: 80-100 recommended
       Set this to be negative if you don't want subsetting
SUBP:  Periods greater than this are fit with the full data,
       periods less than this with the subset data. 
FULL:  1 = Refine displayed estimates with full data.            
       0 = Do not refine with full data before displaying 
MINP:  Minimum period considered:                              
       Halving this quantity almost doubles run-time 
MAXP:  Maximum period considered:                              
       If <=0, maximum period is length of the data
RATE:  Sampling rate in the initial estimation grid: (6 or more) 
       Doubling this quantity almost doubles run-time  
NMIN:  Number of intermediate estimates: 10-20 recommended
FWID:  Width of final grids in standard units(1/span): 1-1.5 recommended  
         Change at your own risk!
FNUM:  Number of values in each final grid: 100-200 recommended  
         Change at your own risk!
SHOW:  Number of period estimates displayed in output

*/

/* --- prototypes ---------------------------------------------------------- */

//#include <stdlib.h>
//#include <math.h>
#include <string.h>
#include <math.h>
#include <libclif.H>
#include "phase.H"

extern "C" {
   int	supdbl_(int*,double*
		,double*,double*,int*,double*,double*,double*,double*);
   int	dqrdc_(double*,int*,int*,int*,double*,double*,double*,int*);
   int	dqrsl_(double*,int*,int*,int*,double*
	       ,double*,double*,double*,double*,double*,double*,int*,int*);
   int	percur_(int*,int*,double*,double*,double*,int*,double*
		,int*,int*,double*,double*,double*,double*,int*,int*,int*);
   int	splev_(double*,int*,double*,int*,double*,double*,int*,int*);
}

void    sort(double *x,int *ii,int n);

/*--- internal functions ---*/
int _ppar( FILE *fp, char *name, double par ) {
   fprintf( fp,"#%-40s: %.4g\n", name, par );
   return 0;
}
int _ppar( FILE *fp, char *name, float par ) {
   fprintf( fp,"#%-40s: %.4g\n", name, par );
   return 0;
}
int _ppar( FILE *fp, char *name, int par ) {
   fprintf( fp,"#%-40s: %d\n", name, par );
   return 0;
}
   
/* --- ReimannPhaser methods ----------------------------------------------- */

void	ReimannPhaser::init( void )
{
   // initialize global data vectors
   times = y = wei = fit = 0;

   // initialize fitting methods' persistent state
   rank = iwrk = 0;
   newt = newy = neww = newf = work = knots = coef = 0;
   res = qraux = dummyp = dummyn = blankp = wy = wxvec = 0;
   wx = 0;

   // initialize memory used locally by phase()
   ord = oo = cc = order = seen = cross = 0;
   ff = mm = FF = tt = TT = yy = YY = ii = rss = rsa = ll = LL = 0;
   ss = SS = initfreq = initmeas = finalfreq = finalmeas = stack = 0;
}

ReimannPhaser::ReimannPhaser( void )
{
   init();
}

ReimannPhaser::~ReimannPhaser( void )
{
   /* 
      this occasionaly causes problems; since i'm only using the phaser 
      in single-shot mode, just don't free anything.  easier than debugging
   */

   //free();
}

ReimannPhaser::ReimannPhaser( const ClifParameters &par )
{
   init();
   load( par );
}

void	ReimannPhaser::free( void )
{
   // release global data vectors
   if ( times ) delete [] times;
   if ( y ) delete [] y;
   if ( wei ) delete [] wei;
   if ( fit ) delete [] fit;
   // release fitting methods' persistent state
   if ( rank ) delete [] rank;
   if ( newt ) delete [] newt;
   if ( newy ) delete [] newy;
   if ( neww ) delete [] neww;
   if ( newf ) delete [] newf;
   if ( work ) delete [] work;
   if ( knots ) delete [] knots;
   if ( coef ) delete [] coef;
   if ( res ) delete [] res;
   if ( qraux ) delete [] qraux;
   if ( dummyp ) delete [] dummyp;
   if ( dummyn ) delete [] dummyn;
   if ( blankp ) delete [] blankp;
   if ( wy ) delete [] wy;
   if ( wxvec ) delete [] wxvec;
   if ( wx ) delete [] wx;

   // release memory used locally by phase()
   if ( ord ) delete [] ord;
   if ( ff ) delete [] ff;
   if ( mm ) delete [] mm;
   if ( FF ) delete [] FF;
   if ( oo ) delete [] oo;
   if ( tt ) delete [] tt;
   if ( TT ) delete [] TT;
   if ( yy ) delete [] yy;
   if ( YY ) delete [] YY;
   if ( ii ) delete [] ii;
   if ( rss ) delete [] rss;
   if ( rsa ) delete [] rsa;
   if ( ll ) delete [] ll;
   if ( LL ) delete [] LL;
   if ( ss ) delete [] ss;
   if ( SS ) delete [] SS;
   if ( cc ) delete [] cc;
   if ( initfreq ) delete [] initfreq;
   if ( initmeas ) delete [] initmeas;
   if ( order ) delete [] order;
   if ( seen ) delete [] seen;
   if ( finalfreq ) delete [] finalfreq;
   if ( finalmeas ) delete [] finalmeas;
   if ( cross ) delete [] cross;

   if ( stack ) delete [] stack;
}

void	ReimannPhaser::load( const ClifParameters &par )
{
   TYPE = par.get( "ReimannPhaser.TYPE", 2 );
   HARM = par.get( "ReimannPhaser.HARM", 1 );
   BASS = par.get( "ReimannPhaser.BASS", 0.0 );
   KNOT = par.get( "ReimannPhaser.KNOT", 20 );
   SUBN = par.get( "ReimannPhaser.SUBN", -100 );
   SUBP = par.get( "ReimannPhaser.SUBP", 4.0 );
   
   //FULL = par.get( "ReimannPhaser.FULL", 0 );
   FULL = par.get( "ReimannPhaser.FULL", 1 );

   //MINP = par.get( "ReimannPhaser.MINP", 0.1 );
   MINP = par.get( "ReimannPhaser.MINP", 0.2 );

   //MAXP = par.get( "ReimannPhaser.MAXP", 100 );
   MAXP = par.get( "ReimannPhaser.MAXP", -1 );

   RATE = par.get( "ReimannPhaser.RATE", 6.0 );
   NMIN = par.get( "ReimannPhaser.NMIN", 20 );
   FWID = par.get( "ReimannPhaser.FWID", 1.0 );
   FNUM = par.get( "ReimannPhaser.FNUM", 100 );
   SHOW = par.get( "ReimannPhaser.SHOW", 15 );
}

int	ReimannPhaser::phase( const CommonLc &lc )
{
  int     i,j,k,l,m;
  double  (ReimannPhaser::*fitcurve)(double);	// pointer to fitting function

/*
   Allocate Space for Data
*/

  free();
  init();

  times = new double [lc.nobs];
  y     = new double [lc.nobs];
  wei   = new double [lc.nobs];
  fit   = new double [lc.nobs];
  ord   = new int [lc.nobs];

  ff    = new double [NMIN];
  mm    = new double [NMIN];
  FF    = new double [NMIN];
  oo    = new int [NMIN];

  tt    = new double [SHOW];
  TT    = new double [SHOW];
  yy    = new double [SHOW];
  YY    = new double [SHOW];
  ii    = new double [SHOW];
  rss   = new double [SHOW];
  rsa   = new double [SHOW];
  ll    = new double [SHOW];
  LL    = new double [SHOW];
  ss    = new double [SHOW];
  SS    = new double [SHOW];
  cc    = new int [SHOW];


  int F, T, S;
  if( sscanf( lc.title, "%*s %*s %d.%d.%d", &F, &T, &S ) != 3 )
     strcpy( idstr, lc.title );
  else
     sprintf( idstr, "%d.%d.%d", F, T, S );

  n = 0;
  for ( i = 0; i < lc.nobs; ++i )
     if ( lc.error[i] > 0 ) {
	times[n] = lc.date[i];
	y[n] = lc.value[i];
	wei[n] = 1 / lc.error[i];
	++n;
     }
  subset = 0;

/* 
   Subset Data 
*/

  fullspan = subspan = times[n-1]-times[0];
  nfull = nsub = n;
  substart=0;
  if(SUBN > 0 && SUBN < n) {
    nsub = SUBN;
    for(i=0;i<=n-nsub;i++) {
      tempspan = times[i+SUBN-1]-times[i];
      if(tempspan < subspan) {subspan=tempspan; substart=i;}
    }
  }

/* 
   Allocate Initial Grid 
*/

  freqmax    = 1.0/MINP; 
  freqchange = 1.0/SUBP;
  if(MAXP > 0) freqmin  = 1.0/MAXP; 
  else {
    if(nsub<nfull && freqchange < 1.0/subspan) freqmin = 1.0/subspan; 
    else freqmin = 1.0/fullspan;  
  }

  if(freqmax<freqmin) {
    fprintf( stderr, "phase: Maximum Period < Minimum Period\n" );
    return 1;
  }

  if(nsub==nfull || freqmax<freqchange) /* All Full Data */
  {
    fullrange = freqmax-freqmin;
    nfreqfull = int( RATE*fullspan*fullrange+1 );
    nfreqsub  = 0;
  }
  else if(freqchange<freqmin)  /* All Subset Data */
  {
    subrange  = freqmax-freqmin;
    nfreqsub  = int( RATE*subspan*subrange+1 );
    nfreqfull = 0;
  }
  else  /* Both types of data in initial grid calculations */
  {
    fullrange = freqchange-freqmin;
    subrange  = freqmax-freqchange;
    nfreqfull = int( RATE*fullspan*fullrange+1 );
    nfreqsub  = int( RATE*subspan*subrange+1 );
  }

  nfreq    = nfreqfull+ nfreqsub;
  initfreq = new double [nfreq];
  initmeas = new double [nfreq];
  order    = new int [nfreq];
  seen     = new int [nfreq];

/* 
   Calculate Initial Grid 
*/

  if(TYPE==1) fitcurve = &ReimannPhaser::trig;
  else if(TYPE==2) fitcurve = &ReimannPhaser::super;
  else if(TYPE==3) fitcurve = &ReimannPhaser::spline;
  fulldata();
  (this->*fitcurve)(1); /* allocate memory in fitcurve */

  for(k=0;k<nfreqfull;k++) {
    initfreq[k] = freqmin + k*fullrange/(nfreqfull-1); 
    initmeas[k] = (this->*fitcurve)(initfreq[k]);
  }
  if(nfreqsub>0) {subsetdata(); m = nfreqfull;}
  for(k=0;k<nfreqsub;k++) {
    initfreq[k+m] = freqmax - (nfreqsub-1-k)*subrange/(nfreqsub-1); 
    initmeas[k+m]   = (this->*fitcurve)(initfreq[k+m]);
  }

/* 
   Find Intermediate Minima 
*/

  sort(initmeas,order,nfreq);
  for(i=0;i<nfreq;i++) seen[i] = 0;
  nbuffer = (int) floor(FWID*0.5*RATE); if(nbuffer<1) nbuffer=1;
  i = 0; j = 0;
  while(i<nfreq && j<NMIN) {
    if(seen[order[i]]==0) /* New minimum! */ {
      ff[j] = initfreq[order[i]];
      for(k=1;k<=nbuffer;k++) { /* wipe out region near minimum */
	if(order[i] >= k)         seen[order[i]-k]=1;
	if(order[i] <= nfreq-1-k) seen[order[i]+k]=1;
      }
      j++; i++; }
    else { /* wipe out neighboring points to boring point */
      if(order[i] >= 1)         seen[order[i]-1] = 1;
      if(order[i] <= nfreq-1-1) seen[order[i]+1] = 1;
      i++; }
  }
  nfound = j;

/* 
   Allocate Final Grids 
*/

  finalfreq = new double [FNUM];
  finalmeas = new double [FNUM];
  final     = nfound*FNUM;

/* 
   Calculate Final Grids
*/

  l=0;
  for(j=0;j<nfound;j++){

     if(ff[j] > freqchange) {
	newrange = (double) FWID/subspan;
	subsetdata();
     }
     else {
	newrange = (double) FWID/fullspan;
	fulldata();
     }   

     newmax = ff[j] + newrange/2;
     newmin = ff[j] - newrange/2;
     if(newmin < freqmin) {
	edge  = 1;
	newmin   = freqmin;
	newrange = newmax - newmin; }
     else if(newmax > freqmax) {
	edge  = 1;
	newmax   = freqmax;
	newrange = newmax - newmin; }
     else edge = 0; 

     for(k=0;k<FNUM;k++) {
	finalfreq[k] = newmin + k*newrange/(FNUM-1);
	finalmeas[k] = (this->*fitcurve)(finalfreq[k]);  
     }

     minloc=0;                
     for(k=1;k<FNUM;k++) if(finalmeas[k]<finalmeas[minloc]) minloc=k; 
     if((minloc>0 && minloc<FNUM-1) || edge==1) {  
	ff[l] = finalfreq[minloc];
	mm[l] = finalmeas[minloc];
	l++; }
  } 
  nfound = l;

/*
   Further Refine Estimates if Requested
*/

  if(FULL==1 && nsub!=nfull) {
     /* sort estimates */
     sort(mm,oo,nfound);
     for(j=0;j<nfound;j++) FF[j] = ff[oo[j]]; 
     for(j=0;j<nfound;j++) ff[j] = FF[j]; 
     if(nfound > SHOW) nfound = SHOW;

     /* cycle over estimates to be refined */
     fulldata();
     for(j=0;j<nfound;j++){
	if(ff[j]>freqchange) {
	   final += FNUM;
	   newrange = (double) FWID/fullspan;
	   newmax = ff[j] + newrange/2;
	   newmin = ff[j] - newrange/2;
	   for(k=0;k<FNUM;k++) {
	      finalfreq[k] = newmin + k*newrange/(FNUM-1);
	      finalmeas[k] = (this->*fitcurve)(finalfreq[k]);  
	   }
	   minloc=0;                
	   for(k=1;k<FNUM;k++) if(finalmeas[k]<finalmeas[minloc]) minloc=k; 
	   ff[j] = finalfreq[minloc];
	   mm[j] = finalmeas[minloc];
	}
     }
  }

/*
   Sort Final Estimates
*/

  sort(mm,oo,nfound);
  for(j=0;j<nfound;j++) FF[j] = ff[oo[j]];  /* sort periods by MOF */
  if(nfound > SHOW)     nfound = SHOW;

/* 
   Averages of Maximum and Minimum y values   
*/

  for(i=0;i<n;i++) fit[i] = y[i];
  sort(fit,ord,n);
  ww=0;WW=0;ee=0;EE=0;
  for(i=0;i<5;i++) {
    ww += pow(wei[ord[i]],2);
    WW += pow(wei[ord[n-i-1]],2);
    ee += y[ord[i]]*pow(wei[ord[i]],2);
    EE += y[ord[n-i-1]]*pow(wei[ord[n-i-1]],2); }
  ee /= ww; EE /= WW;
  /* ee is average of smallest y values, EE of largest y values */

/* 
   Minima, Maxima, WSAR, WRSS, Cycles    
*/

  /* allocate memory */
  cross = new int [lc.nobs];
  stack = new double [lc.nobs];

  for(j=0;j<nfound;j++) {

    /* fit model */

    if(FF[j] > freqchange && FULL==0) subsetdata();
    else fulldata();    
    (this->*fitcurve)(FF[j]);

    /* calc measures of fit */

    rss[j]=0; 
    rsa[j]=0; 
    for(i=0;i<n;i++) {
      rss[j] += pow(wei[i]*(y[i]-fit[i]),2)/n; 
      rsa[j] += fabs(wei[i]*(y[i]-fit[i]))/n; 
    }

    /* maximum and minimum fitted values */

    minloc=0;maxloc=0;
    for(i=1;i<n;i++) {
      if(fit[i]>fit[maxloc]) maxloc=i;
      else if(fit[i]<fit[minloc]) minloc=i; }
    tt[j]=times[minloc]; yy[j]=fit[minloc]; 
    TT[j]=times[maxloc]; YY[j]=fit[maxloc];

    /* Count number of cycles/period in fitted curve */

    /* rescale x and fitted values  */
    for(i=0;i<n;i++) {
      stack[i]=fmod((times[i]-TT[j])*FF[j],1);
      if(stack[i]<0) stack[i]+=1;
      fit[i] = (fit[i]-yy[j])/(YY[j]-yy[j]);
    }
    sort(stack,ord,n);
    /* count crossings */
    oldstate=6;k=0;
    for(i=1;i<n;i++) {
      newstate = (int) ceil(fit[ord[i]]*6);
      if(newstate==0) newstate=1;
      if(newstate==7) newstate=6;
      if(newstate<oldstate)
	for(l=0;l<oldstate-newstate;l++) cross[k++]=0;
      else if (newstate>oldstate)
	for(l=0;l<newstate-oldstate;l++) cross[k++]=1;
      oldstate=newstate;
    }
    /* count cycles */
    l=0; cc[j]=0;
    while(l<k-1) {
      if(cross[l]==cross[l+1]) l++;
      else if(cross[l+1]!=cross[l+2]) l+=2;
      else {cc[j]+=1;l+=2;}
    }
    cc[j] = (cc[j]+1)/2;

    /* calculate average of data close to max/min fitted values */

    ww=0;WW=0;ll[j]=0;LL[j]=0;
    /* find location of max and min in folded space */
    for(i=0;i<n;i++) if(ord[i]==minloc){minloc=i;i=n;}
    for(i=0;i<n;i++) if(ord[i]==maxloc){maxloc=i;i=n;}
    /* loop over neighboring 5 points */
    for(i=-2;i<=2;i++) {
      k = minloc+i;
      l = maxloc+i;
      if(k<0) k+=n; if(k>=n) k-=n;
      if(l<0) l+=n; if(l>=n) l-=n;
      ww += pow(wei[ord[k]],2);
      WW += pow(wei[ord[l]],2);
      ll[j] += y[ord[k]]*pow(wei[ord[k]],2);
      LL[j] += y[ord[l]]*pow(wei[ord[l]],2); 
    }
    ll[j] /= ww; LL[j] /= WW;
    ss[j] = sqrt(1/ww); SS[j] = sqrt(1/WW);
  }

/* 
   Measure-of-fit from Interpolation - WRSS
*/

  for(j=0;j<nfound;j++) {

    /* choose data size */
    if(FF[j] > freqchange && FULL==0) subsetdata();
    else fulldata();

    /* stack times and sort */
    for(i=0;i<n;i++) fit[i]=fmod(times[i]*FF[j],1);
    sort(fit,ord,n);
    for(i=0;i<n;i++) fit[i]=fmod(times[i]*FF[j],1);

    /* calculate interpolation and WRSS */
    ii[j]=0;
    inter(ord[n-1],o=ord[0],ord[1],&est,&var,&wrs); ii[j]+=wrs;
    for(i=1;i<n-1;i++) {
      inter(ord[i-1],o=ord[i],ord[i+1],&est,&var,&wrs); ii[j]+=wrs;
    }
    inter(ord[n-2],o=ord[n-1],ord[0],&est,&var,&wrs); ii[j]+=wrs;
    ii[j] /= n;
  }

    return 0;
}
/*
   End of phase
*/

/* 
   toggle: switch between the full data and the subset data
*/

void	ReimannPhaser::subsetdata( void )
{
   if ( ! subset ) {
      n      = nsub;
      times += substart;
      y     += substart;
      wei   += substart;
      fit   += substart;
      subset = 1;
   }
}

void	ReimannPhaser::fulldata( void )
{
   if ( subset ) {
      n      = nfull;
      times -= substart;
      y     -= substart;
      wei   -= substart;
      fit   -= substart;
      subset = 0;
   }
}

/*
   super:  fit a supersmoother curve to periodic data 
*/

double	ReimannPhaser::super(double f)
{
  int i,two=2;
  double sad,span=0.0,alpha=BASS; 

/* Memory Allocations */
  if ( ! rank ) {
     rank = new int [n];
     newt = new double [n];
     newy = new double [n];
     neww = new double [n];
     newf = new double [n];
     work = new double [7*n];
  }

/* Stack Time and Sort Data */
  for(i=0;i<n;i++) newt[i] = fmod(times[i]*f,1);
  sort(newt,rank,n);
  for(i=0;i<n;i++) {
    newy[i] = y[rank[i]];
    neww[i] = wei[rank[i]];
  }

/* Run Supersmoother */
  supdbl_(&n,newt,newy,neww,&two,&span,&alpha,newf,work);
  for(i=0;i<n;i++) fit[rank[i]] = newf[i];

/* Calculate and Return SAD */
  sad = 0;
  for(i=0;i<n;i++) sad += fabs((y[i]-fit[i])*wei[i])/n;
  return sad;
}
/*
  End of super
*/

/*
   trig:  fit cosines to the data at a certain period    
*/

double	ReimannPhaser::trig(double freq)
{
  int i,j,k,job,ifail,p;
  double RSS,f=2*M_PI*freq; 

  /* Memory Allocation */
  p = 1 + 2*HARM;
  if( ! qraux ) {
     qraux  = new double [p];
     dummyp = new double [p];
     blankp = new double [p];
     res    = new double [n];
     dummyn = new double [n];
     wy     = new double [n];
     wx     = new double* [p];
     wxvec  = new double [p*n];
     for(j=0;j<p;j++) wx[j] = wxvec + j*n;
  }

  /* Design Matrix and Weighted Response */
  for(i=0;i<n;i++) {
    wy[i] = y[i]*wei[i];
    wx[0][i] = wei[i];
    for(k=1;k<=HARM;k++) {
      wx[2*k-1][i] = cos(f*times[i]*k)*wei[i];
      wx[2*k][i] = sin(f*times[i]*k)*wei[i];
    }
  }

  /* Regression */
  ifail = 0;
  dqrdc_(wxvec,&n,&n,&p,qraux,dummyn,dummyp,&ifail);
  ifail = 0;
  job = 111;
  dqrsl_(wxvec,&n,&n,&p,qraux,wy,dummyn,dummyn,blankp,res,fit,&job,&ifail);

  /* Calculate and Return RSS */
  RSS = 0;
  for(i=0;i<n;i++) {
    RSS += res[i]*res[i];
    fit[i] /= wei[i];
  }
  return RSS/n;
} /* End of trig */

/*
      inter: calculate interpolated estimate from neighboring points
*/

void	ReimannPhaser::inter(int i1, int in, int i2
			     , double *EST, double *VAR, double *WRS) 
/* i1 index of left point, in of new point, i2 of right point */
{
  double ca,cb,x1,x2,xn;

  x1=fit[i1]; x2=fit[i2]; xn=fit[in];
  if(x1 > xn) x1 -= 1; /* wrap-around */
  if(x2 < xn) x2 += 1; /* wrap-around */
  ca = (x2-xn)/(x2-x1);
  cb = (xn-x1)/(x2-x1);
  *EST = ca*y[i1] + cb*y[i2];
  *VAR = pow(ca/wei[i1],2) + pow(cb/wei[i2],2);
  *WRS = pow(y[in]-*EST,2);
  *WRS /= pow(1/wei[in],2) + *VAR;
}
/* End of inter */

/*
   spline: fit a periodic cubic spline to data                  
*/

double	ReimannPhaser::spline(double f)
{
  int    i,iopt=-1,m,k=3,nest=KNOT,nknots=KNOT,ier,lwrk;
  double s,wrss; 

/* Memory Allocations */
  m=n+1;
  lwrk=(int)(m*4+nest*23);
  if( ! rank ) {
    rank = new int [m];
    newt = new double [m];
    newy = new double [m];
    neww = new double [m];
    newf = new double [m];
    knots = new double [KNOT];
    coef  = new double [KNOT];
    iwrk  = new int [KNOT];
    work = new double [lwrk];
  }

/* Stack Time and Sort Data */
  for(i=0;i<n;i++) newt[i] = fmod(times[i]*f,1);
  sort(newt,rank,n);
  for(i=0;i<n;i++) {
    newy[i] = y[rank[i]];
    neww[i] = wei[rank[i]];
  }
  newt[n]=newt[0]+1;

/* Initialize knots */
  for(i=1;i<nknots-6;i++)
    knots[i+k]=newt[(int)((double)i/(nknots-6.0)*(n-1.0))];

/* Run Cubic Spline */
  percur_(&iopt,&m,newt,newy,neww,&k,&s,&nest,&nknots,knots,coef,&wrss,work,
	  &lwrk,iwrk,&ier);
  if(ier!=0) {
     fprintf(stderr,"Cubic Spline: ier=%d,freq=%f\n",ier,f);
     return 1;
  }

/* Return Results */
  splev_(knots,&nknots,coef,&k,newt,newf,&m,&ier);
  for(i=0;i<n;i++) fit[rank[i]] = newf[i];
  return wrss/n;
}
/*
  End of spline
*/

int	ReimannPhaser::print( FILE *fp, int verbose ) {

   if( verbose & 256 ) {		// print header only
      _ppar( fp, "ReimannPhaser.TYPE", TYPE );
      _ppar( fp, "ReimannPhaser.HARM", HARM );
      _ppar( fp, "ReimannPhaser.BASS", BASS );
      _ppar( fp, "ReimannPhaser.KNOT", KNOT );
      _ppar( fp, "ReimannPhaser.SUBN", SUBN );
      _ppar( fp, "ReimannPhaser.SUBP", SUBP );
      _ppar( fp, "ReimannPhaser.FULL", FULL );
      _ppar( fp, "ReimannPhaser.MINP", MINP );
      _ppar( fp, "ReimannPhaser.MAXP", MAXP );
      _ppar( fp, "ReimannPhaser.RATE", RATE );
      _ppar( fp, "ReimannPhaser.NMIN", NMIN );
      _ppar( fp, "ReimannPhaser.FWID", FWID );
      _ppar( fp, "ReimannPhaser.FNUM", FNUM );
      _ppar( fp, "ReimannPhaser.SHOW", SHOW );
      return 0;
   }

   if( !fp )
      return SiteErr::errMissing;

   if( !verbose ) {
       fprintf(fp,"%11.7f %8.4f %11.7f %8.4f %11.7f %8.4f %s\n",
	       1/FF[0],mm[0],1/FF[1],mm[1],1/FF[2],mm[2], idstr );
       return SiteErr::errNone;
    }


    /* filename etc. */
    fprintf(fp,"#min_obs max_obs minp maxp rate subp\n");
    fprintf(fp,"%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",ee,EE,
      (double)MINP,1/freqmin,(double)RATE,(double)SUBP);
    fprintf(fp,"#nmin fwid fnum full show nfull nsub fullspan "
	    "subspan nfa_full nf_sub nf_final\n");
    fprintf(fp,"%d %f %d %d %d %d %d %f %f %d %d %d\n",
      (int)NMIN,(double)FWID,(int)FNUM,(int)FULL,
      (int)SHOW,nfull,nsub,fullspan,subspan,nfreqfull,nfreqsub,final);
    fprintf(fp, "#cycles freq     period        "
	    "sup_rsa     sup_rss      int_rss                  "
	    "min_abs       min_fit      min_ave      min_delta     min_se"
	    "                    "
	    "max_abs       max_fit      max_ave      max_delta     max_se\n");
    int j;
    for(j=0;j<nfound;j++) {
       fprintf(fp,"%-2d %12.8f %12.8f %12.8f %12.8f %12.8f             ",
	      cc[j],FF[j],1/FF[j],rsa[j],rss[j],ii[j]);
      fprintf(fp," %12.8f %12.8f %12.8f %12.8f %12.8f               ",
	      tt[j],yy[j],ll[j],yy[j]-ll[j],ss[j]);
      fprintf(fp," %12.8f %12.8f %12.8f %12.8f %12.8f\n",
	      TT[j],YY[j],LL[j],LL[j]-YY[j],SS[j]);
    }


}

/* --- support functions --------------------------------------------------- */

/*
   sort: Phil Spector's sorting routine in C 

     x is a pointer to the array to be sorted
    ii is a pointer to a parallel array of integers to 
       be sorted.  The initialization code in the 
       next line sets ii to 0 ... n-1 - this may or may 
       not be appropriate for your purposes                            
     n is the number of elements in x 
*/

void sort(double *x,int *ii,int n)
{
  int i,ip,iy,lv[16],iv[16],iup,lp;
  register double y;

  for(i=0;i<n;i++)ii[i] = i;
  lv[0] = 0;
  iv[0] = n - 1;
  ip = 0;
  while(ip >= 0)
    {if((iv[ip] - lv[ip]) < 1)
       {ip--; 
	continue; }
     lp = lv[ip] - 1;
     iup = iv[ip];
     y = x[iup];
     iy = ii[iup];
     for(;;)
       {if((iup - lp) < 2)break;
        if(x[++lp] < y)continue;
	x[iup] = x[lp]; 
	ii[iup] = ii[lp];
        for(;;)
           {if((iup-- - lp) < 2)break;
            if(x[iup] >= y)continue;
            x[lp] = x[iup]; 
	    ii[lp] = ii[iup];
	    break; }
        }
     x[iup] = y;
     ii[iup] = iy;
     if((iup - lv[ip]) < (iv[ip] - iup))
       {lv[ip + 1] = lv[ip];
	iv[ip + 1] = iup - 1;
	lv[ip]     = iup + 1; }
     else
       {lv[ip + 1] = iup + 1;
	iv[ip + 1] = iv[ip];
	iv[ip]     = iup - 1; }
     ip++;
    }
} /* End of sort */


int foldLc( CommonLc &lc, CommonLc &plc,
             int *sindex, double *pbuf, double per, double p0 ) {

   /*-------------------------------------------------------------------------*
     foldLc copies lc into plc after converting the dates to phase
     via per and p0 (period and phase0) and resorting.  The header
     info is also copied if the titles and filters are not identical.
     sindex is left with the lc indices sorted in phase order.
    *-------------------------------------------------------------------------*/

   if( lc.nobs <= 0 )
      return SiteErr::errMissing;
 
   if( per <= 0 )
      return SiteErr::errParams;

   /*--- don't recopy header if same name and filter ---*/
   //if( strcmp( lc.title, plc.title ) ||
   //lc.filter != plc.filter )
      lc_hdr_cpy( lc, plc );
 
   plc.resize( lc.nobs );
 
   plc.refJD = lc.refJD;
   if( isnan( p0 ) ) {
      p0 = 0;
   } else {
      plc.refJD -= p0;                  // refJD = phase zero
   }
 
   double *d, *de, *e, *v;
   double *pp = pbuf;
   int n, i, *id;
   d = lc.date;
   double f = 1.0/per;
   for( n = 0; n < lc.nobs; ++n, ++pp, ++d )
      *pp = fmod((*d-p0)*f,1.0);
   
   isort( sindex, pbuf, lc.nobs );
 
   d = plc.date;
   v = plc.value;
   e = plc.error;
   id = plc.obsid;
      
   for( n = 0; n < lc.nobs; ++n, ++d, ++v, ++e, ++id ) {
      i = sindex[n];
      *d = pbuf[i];
      *v = lc.value[i];
      *e = lc.error[i];
      *id = lc.obsid[i];
   }
   return SiteErr::errNone;
}
   
double supersmoothLc( CommonLc &lc, CommonLc &sslc, double alpha ) {

   /*-------------------------------------------------------------------------*
     supersmoothLc 
    *-------------------------------------------------------------------------*/

   static double xv[4096];
   static double yv[4096];
   static double wv[4096];
   static double fv[4096];
   static double work[7*4096];
   double sad,span=0.0;
   int i,two=2;

   if( lc.nobs <= 0 )
      return SiteErr::errMissing;
 
   /*--- don't recopy header if same name and filter ---*/
   //if( strcmp( lc.title, sslc.title ) ||
   //  lc.filter != sslc.filter )
      lc_hdr_cpy( lc, sslc );
 
   sslc.resize( lc.nobs );

   for( i = 0; i < lc.nobs; ++i ) {
      xv[i] = lc.date[i];
      yv[i] = lc.value[i];
      wv[i] = 1.0 / lc.error[i];
   }

  /* Run Supersmoother */
   supdbl_(&lc.nobs, xv, yv, wv, &two, &span, &alpha, fv, work);
   
   sad = 0;
   for( i = 0; i < lc.nobs; ++i ) {
      sad += fabs((yv[i]-fv[i])*wv[i]);
      sslc.date[i] = lc.date[i];
      sslc.error[i] = lc.error[i];
      sslc.value[i] = fv[i];
   }
   printf("%.5f\n", sad/lc.nobs );
   return sad/lc.nobs;
}

int lc_hdr_cpy( const CommonLc &lc1, CommonLc &lc2 ) {
   if( lc2.title )
      delete [] lc2.title;
   lc2.title = new char[strlen(lc1.title)+1];
   strcpy( lc2.title, lc1.title );
   lc2.ra = lc1.ra;
   lc2.dec = lc1.dec;
   lc2.refJD = lc1.refJD;
   lc2.mag0 = lc1.mag0;
   lc2.norm1 = lc1.norm1;
   lc2.norm0 = lc1.norm0;
   lc2.units = lc1.units;
   lc2.filter = lc1.filter;
   lc2.phtpkg = lc1.phtpkg;
   lc2.id = lc1.id;
   return 0;
}
