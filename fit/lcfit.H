#ifndef _lcfit_h
#define _lcfit_h

#include <libsite.H>
#include "libclif.H"

/*--- matrix handling functions ---*/
int matrix( char** &m, int nrows, int ncols );
int matrix( short** &m, int nrows, int ncols );
int matrix( int** &m, int nrows, int ncols );
int matrix( float** &m, int nrows, int ncols );
int matrix( double** &m, int nrows, int ncols );
void free_matrix( char** &m );
void free_matrix( short** &m );
void free_matrix( int** &m );
void free_matrix( float** &m );
void free_matrix( double** &m );
void show_dmatrix( double **m, int R, int C, char *name = NULL );
void show_dvector( double *v, int N, char *name = NULL );
int cgaussj(double **a, int n, double **b, int m); // NRinC modified

/* PeakFinder is the first incarnation of PointFilter, this version is not
   macho specific and can operate on an arbitrary number of passbands. */  
/* written by Mark Pratt April 1995 					*/
/* retained primarily for parameter guessing in LcFit classes below 	*/

class PeakFinder {
public:

   struct element {
      double date;		// date of latest point
      double mdate;		// midpoint of filter points
      double duration;		// time difference between first and last pt
      double max;		// normalized maximum value
      double min;		// normalized minimum value
      double dev;		// sum of chi2 - npoints * penalty (all)
      double pdev;		// sum of chi2 - npoints * penalty (pos)
      double ndev;		// sum of chi2 - npoints * penalty (neg)
      double hipoint;		// largest single pt contrib to deviance (pos)
      double lopoint;		// largest single pt contrib to deviance (neg)
      int Fn;			// Fibonacci number of filter
      int npts;			// number of points in element
      int np;			// number of postive points in element
   };
private:
   LcMoment M;
   int maxn;			// maximum allocated number of measurements
   element *fsingle;		// single point filter
   element *f0;			// filter buffer - even
   element *f1;			// filter buffer - odd
   int fsize;			// size of current filter
   int pfsize;			// size of previous filter
public:
   int npoints;			// current number of measurements
   double penalty;		// chi2 penalty per point (clifpar)
   int nfilters;		// number filters in Fibonacci seq. (clifpar)

   /*--- filter elements of interest ie the finished product ---*/
   element hidev;		// largest filter element (pos)
   element hi2dev;		// 2nd largest filter element (pos)
   element rhidev;		// largest robust filter element (pos)
   element rhi2dev;		// 2nd largest robust filter element (pos)
   element lodev;		// largest filter element (neg)
   element lo2dev;		// 2nd largest filter element (neg)
   element rlodev;		// largest robust filter element (neg)
   element rlo2dev;		// 2nd largest robust filter element (neg)
   
   PeakFinder( ClifParameters &cp );
   ~PeakFinder( void );
   int init( int numpoints );
   int eval( LcSet &lcs, double *norm = 0 );
   int print( FILE *fp = stdout );
};

/* LcFit: multipassband lightcurve fitting 				*/
/* written by Mark Pratt July 1995 					*/
/* fitting adapted from _Numerical Recipes in C_ 			*/
class GeneralLcFit {

protected:

   int maxnobs;			// maximum in any *single* Lc
   int maxnpar;			// maximum number of parameters
   int maxnpb;			// maximum number of passbands

   virtual int eval( double *date, double *fvalue, double **dvdp,
		     int nobs, int pb = -1 ) = 0;
   int load( LcSet &lcs, int npar );

   /*--- Numerical Recipes stuff ---*/
   short mfit;			// Number of fit parameters (!Fixed)
   double **alpha;		// [maxnpar][maxnpar]
   double lamda;

   /*--- NR rewrites ---*/
   int GLFmrqmin( void );
   int GLFmrqcof( double **alpha, double *beta );
   int GLFcovsrt( void );

public:

   enum ParamType { Fixed = 0, Fit, Guess };

   char name[256];
   int Tnobs;			// Total number of observations in LcSet
   int npar;			// Number of parameters
   int nfxpar;			// dof = ndata - (npar - nfxpar)
   double refJD;

   LcSet lcs;	// pointers to data holding members of initializing LcSet
   ClifParameters &cp;
   LcMoment M;
      
   /*--- parameters - buffers managed by fit class ---*/
   double **covar;		// [maxnpar][maxnpar] covarience matrix
   double *par;			// [maxnpar] parameter vector
   double *norm;		// points to lensed normalizations in par[] or null
   double *ptmp;		// [maxnpar] guess parameter vector
   int *ptype;			// [maxnpar] ParamType
   char **pname;		// [maxnpar][64] parameter name

   double *fvalue;		// [maxnobs] fit value
   double **dvdp;		// [maxnpar][maxnobs] dfvalue/dpar

   /*--- convergence criteria ---*/
   double maxdchi2;		// absolute chi2 improvement
   double maxfchi2;		// fractional chi2 improvement
   int maxit;			// max number of iterations
   int nchi;			// number of old chi2 history buffer
   int niter;			// iterations in last minimization (fit flag)
   double pkmin;		// minimum amplification for peak point

   /*--- fit statistics - buffers maintained by fit class ---*/
   double *chi2;		// [maxnpb] fit chi2 of passband
   double *pchi2;		// [maxnpb] peak fit chi2 of passband
   double *mxdchi2;		// largest single point chi2 rel. to baseline
   double *mxdcdate;		// [maxnpb] date of mxdchi2
   int *cnstrt;			// [maxnpb] constraints on passband fit
   int *dof;			// [maxnpb] degrees of freedom in passband
   int *pdof;			// [maxnpb] peak degrees of freedom in passband
   double Chi2;			// cumulative fit chi2
   double PChi2;		// cumulative peak fit chi2
   int DOF;			// cumlative degrees of freedom
   int PDOF;			// cumlative peak degrees of freedom

   GeneralLcFit( ClifParameters &CP = clifpar );
   virtual ~GeneralLcFit( void );

   virtual int init( LcSet &s ) = 0;				// reqd
   virtual int init( LcSet &s, PeakFinder &pf ) = 0;		// not reqd
   virtual int init_par( int parno, double parval, int type = Fit );
   /*--------------------------------------------------------------------------
     iopar: Shows current parameter values and types.  In addition if
     a fit has been performed, parameter errors and fit statistics
     will be displayed.  If prompt is set, the user is prompted for
     new parameter values and types.
     
     mrqfit: NR Levenberg-Marquardt fitting front end.
     If verbose is set, diagnostics are displayed for every iteration
     fitting.
     Convergence is determined using maxdchi2 and maxfchi2.  A buffer
     of chi2 history is compared to the current chi2.  If for
     every chi2 in the buffer
     ochi2[] - chi2 < maxdchi && (ochi2[] - chi2)/chi2 < maxfchi2
     iteration is stopped.  In addition, the routine exits if 
     convergence isn't reached by maxit iterations.  All convergence
     parameters are public members of this class.
     
     renormalize: renormalizes data Lc's, baseline parameters and
     covariance matrix to unit baselines if the fit has identified
     baselines by setting norm;

     ifit: makes calls to iopar, mrqfit and renorm
     
   --------------------------------------------------------------------------*/
   virtual int iopar( int prompt = 0, int stats = 1 );	// parameter setting
   virtual int mrqfit( int verbose );	
   virtual int renormalize( void );
   virtual int ifit( int prompt = 0, int verbose = 0, int renorm = 0 );
   /* FitLc: loads lc with fit values for passband pb at each lc.date  */
   virtual int FitLc( CommonLc &lc, int pb = -1 ) {
      return eval( lc.date, lc.value, 0, lc.nobs, pb ); };
};

class GaussianFit : public GeneralLcFit {

public:
   
   GaussianFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~GaussianFit( void ) { };				// nothing
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf );
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1);
};
   
class MicroFit : public GeneralLcFit {

public:
   
   MicroFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~MicroFit( void ) { };
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf );
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1);
};
   
class MicroBlendFit : public GeneralLcFit {

public:
   
   MicroBlendFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~MicroBlendFit( void ) { };
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf );
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1);
};

class MicroBinSrcFit : public GeneralLcFit {

public:
   
   MicroBinSrcFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~MicroBinSrcFit( void ) { };				// nothing
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf  );
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1 );
};

class MicroFinSrcFit : public GeneralLcFit {
public:
 
   MicroFinSrcFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~MicroFinSrcFit( void ) { };                         // nothing
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf  );
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1);
};

class MicroParFit : public GeneralLcFit {
public:
 
   MicroParFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~MicroParFit( void ) { };                         // nothing
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf  );
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1);
};

class MicroBinSrcFit2 : public GeneralLcFit {

public:
   
   MicroBinSrcFit2( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~MicroBinSrcFit2( void ) { };				// nothing
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf  );
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1 );
};

class MicroYukFit : public GeneralLcFit {
public:
 
   MicroYukFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~MicroYukFit( void ) { };                         // nothing
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf  );
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1);
};

class SineFit : public GeneralLcFit {
public:
 
   SineFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~SineFit( void ) { };                         // nothing
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf  ) { return init( s ); };
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1);
};

class CosineFit : public GeneralLcFit {
public:
 
   CosineFit( ClifParameters &CP = clifpar ) : GeneralLcFit( CP ) { };
   ~CosineFit( void ) { };                         // nothing
   int init( LcSet &s );
   int init( LcSet &s, PeakFinder &pf  ) { return init( s ); };
   int eval( double *date, double *fvalue, double **dvdp, int nobs, int pb=-1);
};

//#endif _lcfit_h
#endif
