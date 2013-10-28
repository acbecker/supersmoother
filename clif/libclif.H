/* written by John Doug Reynolds, April 1998 */

#ifndef _libclif_h
#define _libclif_h

#include <iostream.h>
#include <stdio.h>

/* --- Common Lightcurve Interface core ------------------------------------ */

class CommonLc {
private:
   int			maxobs;		// allocated size of data vectors

public:
   enum Units { Magnitude = 1, Linear };
   enum PhtPkg {
      Undefined = 0, Theory, pSodophot, oSodophot, Sodophot, Daophot, Allframe,
      DIPhot
   };

   char*		title;		// dynamically allocated
   double		ra, dec;	// J2000 coordinates in radians

   double		refJD;		// reference Julian date
   double		mag0;		// zero-point magnitude
   double		norm1;		// ADU = norm1 * value + norm0
   double		norm0;

   short		units;		// enum CommonLc::Units
   short		filter;		// enum Passband::Id
   short		phtpkg;		// enum CommonLc::PhtPkg
   short		id;		// user-defined identification

   int			nobs;		// length of data vectors
   double		*date;		// days from refJD
   double		*value;		// photometric value in units
   double		*error;		// photometric error in units
   int			*obsid;		// observation identifier

   CommonLc( void );
   CommonLc( const char *title, int nobs );
   ~CommonLc( void );

   CommonLc&	operator=( const CommonLc &lc) { copy(lc); return *this; }
   void		copy( const CommonLc&, int header_only = 0 );

   void		erase( void );			// make like new
   int		resize( int nobs );		// extend or truncate
   int		convert( Units u );		// change units
   int		renorm( double n1, double n0 );	// val <-- (val-n0)/n1
   int		print( FILE *fi=stdout, int start=0, int count=-1 ) const;
};

class Passband {
public:
   struct Info {
      double		midpoint;	// Angstroms
      double		bandpass;	// Angstroms
      char		*name;		// printable name
      char		*tag;		// ascii version of id
      short		id;		// enum Passband::Id
   };
   enum Id {
      Invalid = 0, R_Std, V_Std, B_Std, R_MACHO, V_MACHO,
      I_CTIO, R_CTIO, V_CTIO, B_CTIO, U_CTIO,
      I_MJUO, R_MJUO, V_MJUO, B_MJUO, U_MJUO,
      I_UTSO, R_UTSO, V_UTSO, B_UTSO, U_UTSO,
      I_WISE, R_WISE, V_WISE, B_WISE, U_WISE,
      I_MSO30, R_MSO30, V_MSO30, B_MSO30, U_MSO30,
      I_MSO74, R_MSO74, V_MSO74, B_MSO74, U_MSO74,
      I_APO, R_APO, V_APO, B_APO, U_APO,
      I_CTIO15, R_CTIO15, V_CTIO15, B_CTIO15, U_CTIO15,
      RJT_VATT 
   };
   static const Info*	info( const char *tag );
   static const Info*	info( short id );
};

class LcSet {
private:
   short		maxrip;		// allocated size of rip array
   short		numrip;		// number of discarded lightcurves
   CommonLc		**rip;		// array of discarded lightcurves

   short		max;		// allocated size of array
   short		num;		// number of lightcurves
   CommonLc		**lc;		// array of lightcurve pointers

public:
   LcSet( void );
   ~LcSet( void );

   int		count( void ) const { return num; }
   CommonLc*	operator[] ( int index ) const;
   CommonLc*	newlc( void );			// recycled or brand new

   int		append( CommonLc *lc );		// return index of added lc
   int		unhook( const CommonLc *lc );	// remove lc from set
   int		recycle( CommonLc *lc );	// remove and recycle lc
   int		recycle( void );		// recycle entire set
};

class ClifParameters {
public:
   struct Pair { char *name, *value; };

private:
   short		max;
   short		num;
   Pair			*item;
   int			at( const char *name ) const;

public:
   char			*filename;	// name of file currently loaded

   ClifParameters( void );
   ~ClifParameters( void );

   void		reset( void );
   int		load( const char *fname = 0 );
   int		add( const char *name, const char *value );
   const Pair*	operator[] ( int index ) const;
   const char*	operator[] ( const char *name ) const;

   int		get( const char *name, int def ) const;
   double	get( const char *name, double def ) const;
};

// all applications are expected to place clifpar.load() in main
extern ClifParameters clifpar;

/* --- Common Lightcurve Interface packages -------------------------------- */

class BadObsList {
private:
   int 	maxn;				// allocated size of arrays
   int lastr;				// last range, used with ordered ok
public:
   char filename[256];			// complete badobs filename
   short filter;			// Passband::Id
   int nrange;				// number of obsid ranges
   int *min;				// [nrange] min obsid in bad range
   int *max;				// [nrange] max obsid in bad range

   BadObsList( void );
   ~BadObsList( void );

   /*--- mostly internal methods ---*/
   int reset( const int N );
   int reorder( void );
   int match( const short fltr, const char *str );

   /*--- external methods ---*/
   int load( const short fltr, const char *filename, const char *lcdir = 0,
	     const char *object = 0 );
   int ok( const int obsid, const int ordered = 0 );  
};
   
class SodFlags {		 /*--- auxilliary data to CommonLc ---*/
public:
   /* notes on calibration flag:
      0 = uncalibrated
      1 = calibrated with b-r color set to zero
      2 = calibrated with b-r color from the template
      3 = calibrated with b-r color from the running mean
      4 = calibrated with b-r color from the observation	*/

   int			nobs;		// current data length
   short		*crwd;		// crowding parameter
   short		*chi2;		// chi-squared parameter
   short		*calib;		// calibration flag
   float		*airmass;	// airmass
   float		*seeing;	// seeing

   SodFlags( void );
   ~SodFlags( void );

   int 		reset( int size );

private:
   int maxnobs;				// current allocation
};

#ifdef _libsite_h
#ifdef _Macho_h

class SodophotCuts {
public:

   char			*active;	// active production databases
   char			*lcdir;		// local lightcurve cache
   class SodCalibDb	calibdb;	// calibrations database object
   class MachoBadObs	badobs;		// list of bad observations
   int			badobsloaded;	// set by loadbadobs()
   
   // parameters
   char			*badobsFile;
   int			doCuts, verbose, publicCuts, heliocentric;
   int			calibrate, typeMod20, type9ok, typeMax, chi2Max;
   int			crwdMax, cosmMax, missMax, SoP;
   double		errMin, errMax, magMin, chi2magMax, seeingMax;
   double		errScale, errExtra2, refJD;

   ~SodophotCuts( void );
   SodophotCuts( const char *lcdir = 0, const char *active = 0 );
   SodophotCuts( const ClifParameters&
		 , const char *lcdir = 0, const char *active = 0 );

   void	load( const ClifParameters& );
   void set( const char *lcdir, const char *active );
   int	loadbadobs( const char *object = 0 );
   int	clean( LcSet&, SodLcMetaDb& );
   int	clean( const SodophotSet&, int starindex
	       , class CommonLc *rlc, class CommonLc *blc
	       , class SodFlags *rsfp = 0, class SodFlags *bsfp = 0 );
   int	clean( LcSet&, const SodophotSet&, int starindex
	       , class SodFlags *rsfp = 0, class SodFlags *bsfp = 0 );
   int	ok( const SodophotSet::StarObs&, double seeing ) const;
};

#endif /* defined _Macho_h */

class DaophotCuts {
public:
   BadObsList		badobs;
   char			*lcdir, *badobsfile;
   int			doCuts, normMin;
   double		sharpMin, sharpMax, aspectMax, quadErrMax, dMagChi2Max;
   double		errScale, errExtra2;
   double		refJD;

   DaophotCuts( const char *lcdir = 0 );
   DaophotCuts( const ClifParameters &params, const char *lcdir = 0 );
   ~DaophotCuts( void );

   void		load( const ClifParameters &params );
   void		set( const char *lcdir );
   int		clean( CommonLc&, const DaophotLc& );
   int		clean( LcSet&, const DaophotSet& );
};

#endif /* defined _libsite_h */

/* commented this out to make fazer subdir compile, its defined in phase.H */

/*
class ReimannPhaser {
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
   double		ee,EE,ww,WW,est,var,wrs,*stack,*ff,*mm,*FF,*tt,*TT;
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

private:
   // persistent state for fitting methods
   int			*rank, *iwrk;
   double		*newt,*newy,*neww,*newf,*work,*knots,*coef;
   double		*res,*qraux,*dummyp,*dummyn,*blankp,*wy,**wx,*wxvec;
   void		init( void );
};
*/

/* --- general utilities --------------------------------------------------- */

char*	duplicate( const char *str );
int	isort( int index[], const double value[], int length );
int	isort( int index[], const int value[], int length );
int	fold( class CommonLc&, const CommonLc&, double period, double date0 );
double	heliocentric( double juliandate, double ra, double dec );

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

/* --- Mark Pratt's Utilities ---------------------------------------------- */

int	LcMerge( const CommonLc &lc1, const CommonLc &lc2, CommonLc &lc
		 , const ClifParameters &par = clifpar );
void	LcMerge( LcSet &set, const ClifParameters &par = clifpar );

class LcMoment {
private:
   int ok;			// successful eval completed
   int sorted;			// successful sorting completed
   short maxn;			// max size of sindex
   short _size;			// current size of sindex
   int nlo, nhi;		// upper and lower index limits
   int *sindex;			// sorted index array
   
public:

   int N;			// Sum 1
   double S;			// Sum 1/err^2
   double Sv;			// Sum value/err^2
   double Sdv2;			// Sum dvalue^2/err^2
   double Sdv3;			// Sum dvalue^3/err^2
   double Sdv4;			// Sum dvalue^4/err^2
   double Se;			// Sum error

   LcMoment( void );
   ~LcMoment( void );
   int reset( int n );
   int eval( CommonLc &lc, double ptilelo = 0.0, double ptilehi = 1.0);
   int pindex( double percentile = 0.5 );
   double moment( int n );	// moments are per d.o.f.
   double average( void );	// Sv/S
   double chisq( void );	// Svv - Sv^2/S	
   double skewness( void );	// moment(3)/moment(2)^1.5
   double kurtosis( void );	// moment(4)/moment(2)^2
   int print( FILE *fp = stdout, int verbose = 0 );
};

/* --- SM Ploting Utilities ------------------------------------------------ */

      /*--- calls with CommonLc will expand passed in limits ---*/
int LcDataLimits( CommonLc &l, float &tlo, float &thi,
                  float &vlo, float &vhi, float &elo, float &ehi );
int LcDataLimits( CommonLc &l, float &tlo, float &thi, float &vlo, float &vhi);

      /*--- calls with LcSet will reset passed in limits ---*/
int LcDataLimits( LcSet &s, float &tlo, float &thi,
                  float &vlo, float &vhi, float &elo, float &ehi );
int LcDataLimits( LcSet &s, float &tlo, float &thi, float &vlo, float &vhi );


class LcPlot {
private:
   int maxn, nobs;
   int maxp;				// maximum number of passbands
   int ok;
   float *xb, *yb, *eb;
   int loadlc( CommonLc &lc );		// used by exclusively by plot
public:
   /*--- overwrite these for custom configuration ---*/
   char devstr[1024];			// sm device
   char plabel[1024];			// plot label
   char xlabel[1024];
   char ylabel[1024];
   float xnot, ynot;			// notation placement [0,1] in window
   float ptype[2];
   float xlimlo, xlimhi;		// global x limits
   float ylimlo, ylimhi;		// global y limits 
   float *ylim[2];			// limits of individual lightcurve

   LcPlot( void );
   ~LcPlot( void );
   int resize( int size );
   int config( LcSet &s );		// use LcSet info to configure 
   int draw( LcSet &s );		// plot only using prev. config
   int print( LcSet &s, const char *dev, const char *filename = 0 );
};

#endif /* ! defined _libclif_h */
