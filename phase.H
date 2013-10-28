#ifndef _phase_h
#define _phase_h

#include <libsite.H>
#include <libclif.H>

/*--- ReimannPhaser: written by James Reimann, August 1994 */
/*--- repackaged by John Doug Reynolds, July 1996 */   

class ReimannPhaser {
private:
   // persistent state for fitting methods
   int			*rank, *iwrk;
   double		*newt,*newy,*neww,*newf,*work,*knots,*coef;
   double		*res,*qraux,*dummyp,*dummyn,*blankp,*wy,**wx,*wxvec;
   void		init( void );

public:
   // global parameters

   int			TYPE,HARM,KNOT,SUBN,FULL,NMIN,FNUM,SHOW;
   double		BASS,SUBP,MINP,MAXP,RATE,FWID;

   // global variables

   char 		idstr[64];	// usually field.tile.seq
   int			n;		// number of observations
   int			subset;		// currently in subset mode
   double		*times;		// vector of observation times
   double		*y;		// vector of signals
   double		*wei;		// vector of regression weights
   double		*fit;		// vector of fitted values

   // variables used locally by phase()

   //			subsetting
   int			substart,nsub,nfull;
   double		tempspan,fullspan,subspan;

   //			initial grid
   int			nfreq,nfreqfull,nfreqsub,final,nfound;
   double		freqmax,freqmin,freqchange,fullrange,subrange;
   double		*initfreq,*initmeas;

   //			final grids
   int			edge,nbuffer,*seen,*order;
   double		newrange,newmax,newmin,*finalfreq,*finalmeas;

   //			final estimates
   int			minloc,maxloc,o,oldstate,newstate,*oo,*ord,*cc,*cross;
   double		ee,EE,ww,WW,est,var,wrs,*ff,*stack,*mm,*FF,*tt,*TT;
   double		*yy,*YY,*ii,*rss,*rsa,*ll,*LL,*ss,*SS;

   // support functions

   double	super(double f);
   double	trig(double f);
   double	spline(double f);
   void		inter(int i1,int in,int i2,double*est,double*var,double*wrs);
   void		subsetdata(void);
   void		fulldata(void);

   // primary user interface routines

		ReimannPhaser( void );
		ReimannPhaser( const ClifParameters& );
		~ReimannPhaser( void );

   void		free( void );
   void		load( const ClifParameters& );
   int		phase( const CommonLc& );
   int 		print( FILE *fp = stdout, int verbose = 0 );
};

double supersmoothLc( CommonLc &lc, CommonLc &sslc, double alpha );
int foldLc( CommonLc &lc, CommonLc &plc,
             int *sindex, double *pbuf, double per, double p0 );

int lc_hdr_cpy(const CommonLc &lc1, CommonLc &lc2);

/* SSTab: tabulated supersmoother 					*/
/* written by Mark Pratt  Oct 1996 					*/
class SSTab {
private:
   int maxnobs;						// current allocation
public:
   int ndata;						// data vector length
   int ok;						// current state
   double F;						// 1 / period
   double *dv, *xv, *yv, *wv;				// data vectors
   double *fv, *dfv;					// smoothed vectors
   double *work;					// ss work space
   double BASS;
   int *pindex;						// phase -> date
   int *dindex;						// date -> phase
   int *nexti;						// data index

   SSTab( void );
   ~SSTab( void );
   int reset( int npts );
   int init( const CommonLc &lc, double period );
   int init( const CommonLc &lc, double period, double p0);
   double ss( double d );				// project ss 
   /* Warning: dss is a derivative w.r.t. phase not date */ 
   int eval( double d, double *ss, double *dss );	// project ss, ss'
   double chi2( void );					// chi2/df
   double fom( void );					// figure of merit/df
};

/* DPhase: periodic lightcurve timing 					*/
/* written by Mark Pratt Oct 1996 					*/
class DPhase : public SSTab {
public:
   int obsperfit;			// observations per phase fit
   double tunits;			// 1 or days --> AU conversion
   double maxfrac;			// maximum fraction of timespan/bin
   CommonLc pLc;			// phased data Lc
   CommonLc ssLc;			// phased supersmoothed Lc
   CommonLc dpLc;			// phase evolution Lc (date)

   DPhase( ClifParameters &CP = clifpar );
   ~DPhase( void );
   int derive( const CommonLc &dataLc, double period );
   int alphabet(  double *alpha, double *beta, double *chi2, double dp,
		  int ilo, int ihi );
   int fitphase( int ilo, int ihi, int ndp ); 
};


class FitPhase : public SSTab {
public:
   double period, dperiod;
   CommonLc jLc, hLc, kLc;
   SSTab    jSs, hSs, kSs;

   FitPhase( );
   ~FitPhase( void );
   int alphabet(  double *alpha, double *beta, double *chi2, double dp );
   int fitphase( CommonLc &jLc, CommonLc &hLc, CommonLc &kLc ); 
   int fitphase( CommonLc &Lc);
};


#endif
