/* written by John Doug Reynolds, April 1998 */

#ifndef _libsite_h
#define _libsite_h

/*         To most fully utilize this header, first include Macho.h          */

#include <sys/types.h>
#include <iostream.h>
#include <stdio.h>
#include <time.h>
#include "prototypes.h"

/* --- Miscellaneous ------------------------------------------------------- */

extern double Infinity, NaN;
extern "C" void initialize_site_globals( void );

// convert between time_t and macho date (days since 2 Jan 1992)
extern double	machotime( int time );
extern int	timemacho( double days );

// convert between time_t and Julian date
extern double	juliantime( int time );
extern int	timejulian( double jd );

// convert from macho RA and DEC to static string
extern char *	macho_ra_string( int ra, char delim = ' ', int prcsn = 3 );
extern char *	macho_dec_string( int dec, char delim = ' ', int prcsn = 2 );

// convert from string to macho RA and DEC
extern int	macho_ra_string( const char *str );
extern int	macho_dec_string( const char *str );

// construct names for things in a local cache (lcdir)
class SodophotSet;
extern char *	lcdir_star_dir( const char *lcdir, const char *starid
				, int mode = 'r', char *buf = 0 );
extern char *	lcdir_sodlc_file( const char *lcdir, const char *starid
				  , int source, char *buf = 0 );
extern int	lcdir_sodlc_find( const char *lcdir, const char *starid
				  , int source, class SodophotSet &ss );

// wrappers for complicated POSIX signals calls
typedef void (*sigf_t)(int);
extern int	signal_catch( int signal, sigf_t handler, int flags = 0 );
extern int	signal_unblock( int signal );
extern int	signal_block( int signal );

/* --- SiteErr ------------------------------------------------------------- */

class SiteErr {
public:

   enum {
      errNone = 0, errProgram, errParams, errEnviron, errMemory,
      errAccess, errFormat, errMissing, errVersion, errMath, errCorrupt,
      MaxErrs = 100
   };

   int		errcode;

   SiteErr( int error = errNone ) { errcode = error; }
   virtual ~SiteErr( void ) {}
   operator int() const { return errcode; }
   SiteErr& operator= ( int error ) { errcode = error; return *this; }
   virtual int perror( const char *label = 0 ) const;
};

/* --- ExplicitTM ---------------------------------------------------------- */

class ExplicitTM : public tm {
private:
   char			*out;
   int			maxout;
   int			format( const char *fmt, int &len );
   int			scan( const char* &str, const char *fmt );

public:
   enum { ZoneMax = 64 };

   char			zone[ZoneMax];
   char			zone3[3];
   int			gmtoff;

   ExplicitTM( void );
   ~ExplicitTM( void );

   int		reset( const char *zone = 0 );
   int		reset( const time_t &clock, const char *zone = 0 );
   const char*	etime( const time_t &clock, const char *zone = 0 );
   const char*	sscan( const char *str, const char *fmt );
   const char*	format( const char *fmt = 0 );
   int		get( time_t &clock );
};

/* --- ChunkToken ---------------------------------------------------------- */

#ifdef _Macho_h

class ChunkToken {
private:
   const static short	red[8][8];	// East Red map
   const static short	blue[8][8];	// East Blue map
   short		tok;		// encodes sop, clr, row, col

   ChunkToken( int t ) { tok = short(t); }

public:
   ChunkToken( void ) { tok = -1; }
   ChunkToken( SideOfPier sop, short id );
   ChunkToken( SideOfPier sop, Color clr, int row, int col );
   ChunkToken( SideOfPier sop, short chip, int row, int col );

   SideOfPier	sop( void ) const { return tok & 0200 ? West : East; }
   Color	clr( void ) const { return tok & 0100 ? Blue : Red; }
   int		row( void ) const { return (tok & 070) >> 3; }
   int		col( void ) const { return tok & 07; }

   short	chip( void ) const;
   short	id( void ) const;

   int		operator! ( void ) const { return tok < 0; }
   ChunkToken	operator- ( void ) const { return ChunkToken(tok ^ 0200); }
   ChunkToken	operator~ ( void ) const { return ChunkToken(tok ^ 0100); }

   ChunkToken	incr( int row, int col ) const;
   ChunkToken	operator+ ( SideOfPier s ) const;
   ChunkToken	operator+ ( Color c ) const;
};

#endif /* _Macho_h */

/* --- SodLcDb ------------------------------------------------------------- */

class SodLcIndex : public SiteErr {
private:
   int			size;
   char			*map;

public:
   enum { MaxTile = 65536 };

   struct LcIndex {
      short	field, chunk;
      u_short	seq1, seq2;
      u_short	index, count;
   };

   SodLcIndex( void ) { map = 0; }
   SodLcIndex( const char *path, const char *name = 0 );
   ~SodLcIndex( void ) { close(); }

   int		open( const char *path, const char *name = 0 );
   int		close( void );

   int		count( int tile ) const;
   LcIndex*	operator[] ( int tile ) const;
};

class SodLcDb : public SiteErr {
public:
   struct Header {
      int		_rec1_start;
      char	comment[80];
      int		_rec1_end;
      int		_rec2_start;
      int	iwjsver;
      int		_rec2_end;
      int		_rec3_start;
      int	chunk;
      int	tileno[10];
      int		_rec3_end;
      int		_rec4_start;
      int	nobs;
      int	nstars;
      int		_rec4_end;
      int		_rec5_start;
      char	comment2[80];
      int		_rec5_end;
      int		_rec6_start;
      float	spare[20];
      int		_rec6_end;
      int		_rec7_start;
      float	spare2[20];
      int		_rec7_end;
   };

   struct Observation {
      int		_rec_start;
      int	obsno;
      int	ipier;
      float	time;
      float	seeing[2];
      float	avsky[2];
      float	airmass;
      float	sodver;
      float	tempver;
      float	tempobs[2];
      float	dxobs;
      float	dyobs;
      int	isodtime;
      float	dbvers;
      int	idbtime;
      int	lot;
      int	checklist;
      int	ffflags;
      float	spare;
      int		_rec_end;
   };

   struct Photometry {
      short	rmag, rerr, rflag1, rflag2, rflag3;
      short	bmag, berr, bflag1, bflag2, bflag3;
   };

   struct Lightcurve {
      int			_rec1_start;
      int		tile;
      int		seq;
      int		nobs1;
      int		nsobs[2];
      int			_rec1_end;
      int			_rec2_start;
      int		temptype[2];
      float		x[2];
      float		y[2];
      int			_rec2_end;
      int			_rec3_start;
      int		ra;
      int		dec;
      float		spare[4];
      int			_rec3_end;
      int			_rec4_start;
      SodLcDb::Photometry
	                phot[1];		// see sizeof_Lightcurve
      int			_rec4_end;
   };

private:
   int			filesize;
   size_t		sizeof_Lightcurve;
   char			*lc0;

public:
   char			*fname;
   Header		*hdr;
   Observation		*obs;

   SodLcDb( void );
   SodLcDb( const char* fname );
   SodLcDb( const char* fname, int nstars, int nobs );
   ~SodLcDb( void ) { close(); }

   int		open( const char* fname );
   int		open( const char* fname, int nstars, int nobs );
   int		close( void );

   Lightcurve*	operator [] ( int index )
      { return (Lightcurve*) (lc0 + index * sizeof_Lightcurve); }

   Lightcurve*	find( int tile, int seq );
   Lightcurve*	find( int tile, int seq, int index, int count );
};

class SodLcMetaDb : public SiteErr {
private:
   struct DbTOC {
      int		obsid;
      int		rtime;
      int		segment;
      int		obsindex;
   };

   int			ndb;		// number of segments
   SodLcDb		*db;		// array of database segments
   int			nobs;		// number of toc entries
   DbTOC		*toc;		// table of contents
   int			tile, seq;	// current star to follow
   int			primary;	// primary database segment
   int			segment;	// current database segment
   int			index;		// current toc index

   static int	compare( const void *i, const void *j );

public:
   // pass-through interface to current database segment
   int			field;
   char			*fname;
   SodLcDb::Header	*hdr;
   SodLcDb::Observation	*obs;
   SodLcDb::Lightcurve*	operator[] ( int starindex ) const;
   SodLcDb::Lightcurve* find( int tile, int seq, int index, int count ) const;
   SodLcDb::Lightcurve* find( int tile, int seq ) const;

   // these are maintained by reset() and next()
   SodLcDb::Lightcurve	*star;		// current lightcurve
   int			obsindex;	// current observation index

   SodLcMetaDb( void );
   ~SodLcMetaDb( void );

   int		open( const char *path, int field, int chunk, int wjstype );
   int		close( void );

   int		obscount( void ) const;		// max number of observations
   int		starcount( void ) const;	// max number of stars
   int		setdb( int segment = -1 );	// change current database

   int		reset( int tile, int seq );	// reset to first obs
   int		reset( int starindex );		// reset to first obs
   int		next( void );			// proceed to next obs
};

/* --- DaophotLc ----------------------------------------------------------- */

class DaophotLc : public SiteErr {
public:
   struct Photometry {
      int			date;
      int			obsid;
      float			mag, err, sky, xpix, ypix, airmass;
      float			xFWHM, yFWHM, chi2, sharpness;
      float			dMag, dMagErr, dMagChi2;
      short			normstars;
      short			exposure;
   };

   int				nobs;
   char				*object;
   char				*filter;
   DaophotLc::Photometry	*phot;

   DaophotLc( void );
   ~DaophotLc( void ) { erase(); }

   int		extract( const char *path, const char *object
			 , const char *filter, int verbose = 0 );
   int		read( const char *lcdir
		      , const char *object, const char *filter );

   void		erase( void );
   int		print( FILE *fi = stdout
		       , int skip = 0, int count = -1, int flags = 3 ) const;
   int		write( const char *lcdir );
   int		write( int fd );
   int		read( int fd );
};

class DaophotSet {
private:
   short		max;		// allocated size of array
   short		num;		// number of lightcurves
   DaophotLc		**lc;		// array of lightcurve pointers

public:
   DaophotSet( void );
   ~DaophotSet( void );

   int		count( void ) const { return num; }
   int		append( DaophotLc *lc );
   int		unhook( const DaophotLc *lc );
   DaophotLc*	operator[] ( int index ) const;

   void		erase( void );
   int		extract(const char *path, const char *object, int verbose = 0);
   int		read( const char *lcdir, const char *object );
   int		write( const char *lcdir );
};

/* --- LcStat -------------------------------------------------------------- */

class LcStat : public SiteErr {
private:
   char 	*map;
   int		_tile;
   void		initialize( void );

public:
   struct MicroFilter  {		// microlensing level 1 stats 
      float	filmax;
      float	filtime;
      int   	nhipoints;
   };
   
   struct MonoChromePhot {
      float	mag;		// Median of good measurement: rmagave
      float	errave;		// average error?
      int	N;		// number of good points? 
      float	sig;	   
      float	chi2w;
      float	qw;
      int	Nr;		// number of robust measurements?
      float	erraver;	// robust average error?
      float	chi2r;
      float	qr;
      int	s5;
      int	s5d;
      float	magaveu;	// median of all measurements: rmagaveu 
      float	autoc;		// Welch: auto correlation? 
      float 	win[10];	// Welch: correlation windows 
      float	decorcon; 	// decorrelation against seeing 
      float	decorslope; 	// decorrelation against seeing 
      /*-- auxillary sodophot measure info (normally in sodstat.out1 file) --*/
      float	dskyave;
      float	dskysig;
      float	psfave;
      float	psfsig;
      float	crdave;
      float	cosmave;
      float	missave;
      MicroFilter  m1stat[3];		// statdb.h: #define NUMFILTERS 3
   };
   
   struct ColorPhot {	   	  	// two color info */
      int	N;			// number of simult red/blue points
      float	sig;
      float	chi2w;
      float 	ave;
      float	probspear;		// rank red/blue correlation coeff.
      float	cross;			// Welch: rb cross correlation 
   };

   struct Rec {
      ushort 	tile;
      ushort	seq;
      int	RA;
      int 	Dec;
      int	cat;
      int	KGstatnum;
      int	KGvarnum;
      int	npoints;
      int	aversion;
      int 	lot;
      MonoChromePhot 	red;
      MonoChromePhot 	blue;
      ColorPhot 	br;
   };

   int 		size;
   char		*fname;
   int		Nseq;

   LcStat( void ) { initialize( ); };			// zero pointers
   ~LcStat( void ) { clear( ); };			// deallocate
   void clear( void );					// deallocate & zero
   int open( char *path, int tile );			// open for reading
   int open( char *path, int tile, int nrec );		// open for read/write
   Rec * operator[] (int sequence) const;
   static int print( Rec &rec, FILE *fp = stdout );
   int print( int seq, FILE *fp = stdout );

};

/* --- SodophotSet --------------------------------------------------------- */

// this is copied from MachoBB.h because that header is not portable
#define		FlatFielded	0x000010

struct ExtractLcOpts {
   ExtractLcOpts( void );		// assigns default settings
   const char*		active;		// path to active databases
   const char*		cache;		// path to database cache
   const char*		vmfile;		// optional output file name
   int			verbose;	// print diagnostic messages
   int			tolerant;	// skip inaccessible databases
   int			minobs;		// reject earlier observations
   int			maxobs;		// reject later observations
   int			quota;		// sufficient number of obs
   int			merge;		// limit on database merges
   int			unzip;		// limit on uncompressions
   int			fetch;		// limit on nearline fetches
   int			field;		// mandatory field number
   int			tile;		// mandatory tile number
};

class SodophotSet : public SiteErr {
public:
   struct Coordinates {
      float			x, y;		// chip coordinates
   };
   struct Dimensions {
      int			s, o, c;	// stars, observations, chunks
   };
   struct Header {
      short			field;		// field number
      unsigned short		tile;		// tile number
      struct Dimensions		count;		// occupied dimensions
      int			source;		// source database code
      int			incept;		// template inception time_t
   };
   struct Star {
      unsigned short		seqn;		// star sequence number
      unsigned char		wcid[2];	// canonical chunkid [color]
      struct Coordinates	pix[2];		// template coords [color]
      short			mag[2];		// template milli-mag [color]
      int			ra, dec;	// milli-arcseconds
   };
   struct Observation {
      int			date;		// end of exposure time_t
      int			obsid;		// observation number
      int			checklist;	// pipeline status flags
      short			exposure;	// seconds of exposure
      short			pierside;	// East or West
      float			airmass;	// airmass
   };
   struct Photometry {
      short			mag, err;	// milli sod magnitude
      short			flag[3];	// packed photometry flags
   };
   struct ChunkObs {
      int			tobsid;		// template observation number
      signed char		ampid;		// amplifier id number
      signed char		ccdflips;	// sky orientation flags
      unsigned char		offset[2];	// amp offset in FITS [x,y]
      float			matrix[6];	// coordinate transform matrix
      float			seeing;		// sqrt( PSF area )
      float			avesky;		// average sky level
   };
   struct StarObs {
      float			mag, err;	// sod magnitude units
      short			dsky, type;	// unpacked flag 1
      short			crwd, chi2;	// unpacked flag 2
      short			mpix, cosm;	// unpacked flag 3
      short			ampid;		// corrected amplifier id
      struct Coordinates	pix;		// approximate coordinates
   };

protected:
   int				bfv;		// binary format version
   char*			map;		// memory map
   struct Dimensions		max;		// allocated dimensions
   struct Photometry*		phot;		// phot[ nseq ][ nobs ][ 2 ]

   int		mapsize( void );

public:
   struct Header*		hdr;		// hdr[ 1 ]
   struct Star*			star;		// star[ nseq ]
   struct Observation*		obs;		// obs[ nobs ]
   struct ChunkObs*		chnk[128];	// chnk[ 128 ][ nobs ]

		SodophotSet( void );
		~SodophotSet( void );

   void		close( void );
   int		allocate( int nseq, int nobs );
   int		allocate( int cid );

   int		open( const char *file, char mode = 'r' );
   int		open( const char *file, int nseq, int nobs, int nchk );
   int		write( const char *file, int count, const int *sindex ) const;
   int		write( const char *file ) const;

   int		extract( SodLcMetaDb& );
   int		extract( ExtractLcOpts, int count, const ushort *seqn );
   int		extract( ExtractLcOpts, int skip=0, int count=0, int inc=1 );

   void		pack( const StarObs&, int sindex, int oindex, int color );

   int		version( void ) const;
   int		find( unsigned short seqn ) const;
   char*	starid( int starindex, char *buf = 0 ) const;

   // character codes for coords: 'c'=chip, 'a'=amp, 'f'=FITS, 0=none
   void		unpack( StarObs&, int starindex
			, int obsindex, int color, char coords = 0 ) const;

   int		printlc( FILE*, int starindex, char coords = 'a'
			 , int obsindex = 0, int count = -1 ) const;
   int		printso( FILE*, int starindex, int obsindex, char coords = 'a'
			 , int rflag = -1, int bflag = -1 ) const;

   struct Dimensions	maximum( void ) const;
   struct Photometry	( *operator[] ( int starindex ) const ) [2];
};

int merge( SodophotSet &sum, const SodophotSet &s1, const SodophotSet &s2 );
int sodset_byte_reverse( const char *dstfile
			 , const char *srcfile, int verbose = 1 );
int excerpt( SodophotSet &dst, const SodophotSet &src
	     , int scount, const int *sindex, int ocount, const int *oindex );

/* --- SodCalibDb ---------------------------------------------------------- */

#ifdef _Macho_h

class SodCalibrator {
public:
   int			field;		// negative for east-style template
   int			rwcid;		// red west chunk id
   int			obsid;		// template observation number
   int			npsf;		// number of psf stars

   float		airmass;	// template airmass
   float		a0, a1;		// V transform: zero pt, color coef
   float		b0, b1;		// R transform: zero pt, color coef
   float		co;		// chunk offset
   float		mc;		// psf mean color

   float		bje, bjw, bjo;	// blue jitter coefficients
   float		V0, V1;		// V transform: linear coefficients
   float		R0, R1;		// R transform: linear coefficients

   SodCalibrator( void );

   /* notes on conversion functions:
      mag:	raw Sod instrumental magnitude
      sop:	pier-side (East or West) of mag
      bmr:	B-R color appropriate for mag		*/

   // Sod instrumental blue magnitude with jitter removed
   float B( SideOfPier sop, float bmag, float bmr ) const
   {		return bmag + (sop == East ? bje : bjw) * (bmr - bjo); }

   // Standard V magnitude
   float V( SideOfPier sop, float bmag, float bmr ) const
   {		return B(sop,bmag,bmr) + V0 + V1 * bmr; }

   // Standard R magnitude
   float R( float rmag, float bmr ) const
   {		return rmag + R0 + R1 * bmr; }
};

class SodCalibDb {
public:
   struct Chip {
      int		obsid;		// template observation number
      float		airmass;	// template airmass
      float		a0, a1;		// V transform: zero pt, color coef
      float		b0, b1;		// R transform: zero pt, color coef
   };
   struct Chunk {
      int		rwcid;		// red west chunk id
      int		npsf;		// number of psf stars
      float		co;		// chunk offset
      float		bj;		// blue jitter coefficient
      float		mc;		// psf mean color
   };

   int			field;		// negative for east-style template
   struct Chip		chip[4];	// index by red chip number
   struct Chunk		chunk[64];	// index by (rwcid + 8) & 63

   SodCalibDb( void );

   /* notes:
      The integer returned by these functions is a SiteErr.
      The second read() constructs the filename from the given field and
      looks first in lcdir, if not null, and then in active, if not null */

   int		read( const char *filename );
   int		read( int field, const char *lcdir, const char *active );
   int		get( class SodCalibrator&, int rwcid ) const;
   int		print( FILE *fi = stdout ) const;
};

#endif /* _Macho_h */

/* --- MachoBadObs --------------------------------------------------------- */

class MachoBadObs {
public:
   struct Range {
      int		min, max;	// inclusive observation range
      int		amps;		// each bad amp has its bit set
   };

   struct Range		*list;		// ordered, non-overlapping ranges
   int			count;		// number of ranges in use

		MachoBadObs();
		~MachoBadObs();

   /* notes:

      The integer returned by load() is a SiteErr.  find() returns the
      index of the range containing obsid, or -1.  Both forms of amps()
      return the amps value associated with obsid.  The first form uses
      index as the starting point for a sequential search, and leaves
      it suitably initialized for the next call.  index need not be
      initialized before the first call, though it does make sense to.
      Invoke normalize() after any manual alteration of the ranges.	*/

   void		normalize( void );
   void		extend( int capacity );
   int		load( const char *filename, int verbose = 1 );
   int		find( int obsid ) const;
   int		amps( int obsid, int &index ) const;
   int		amps( int obsid ) const;
   int		capacity( void ) const;

private:
   int			maxcap;		// allocated length of list
};

/* --- StarNote ------------------------------------------------------------ */

class StarNote {
private:
   int			buflen;		// allocated length of buffer

public:
   short		field;		// star field identifier
   unsigned short	tile;		// star tile identifier
   unsigned short	seqn;		// star sequence identifier
   unsigned short	length;		// occupied length of buffer
   char*		buffer;		// memory for arbitrary data
   char			label[31];	// data format identifier
   char			user[9];	// login name of creator
   int			date;		// creation time stamp

		StarNote( void );
		~StarNote( void );
		StarNote& operator=( const StarNote &note );

   int		set_length( int length );	// change note length safely
   int		set_label( const char *label );	// change format label safely
   int		set_user( void );		// identify process owner

   int		print( FILE *fi = stdout ) const;
   int		scan( FILE *fi );
};

extern const char *category_source_menu[];
extern const char *category_name_menu[];
extern const char *startag_alert_menu[];

extern const char*
select_from_menu( const char *label, const char *list[], FILE *fi = stdout );

extern int
construct_category_note( StarNote &note, short field, ushort tile
			 , ushort seqn, const char *src, const char *cat );
extern int
construct_startag_note( StarNote &note, short field, ushort tile
			, ushort seqn, const char *alert, int cameos );

/* --- FitsHeader ---------------------------------------------------------- */

class FitsHeader : public SiteErr {
private:
   char			(*card)[80];	// array of cards
   int			pages;		// number of pages allocated
   int			end;		// index of END card

   int		grow( int pages );

public:
   char			buffer[71];	// used by formatting methods
   int			cursor;		// index of current active card

		FitsHeader( int pages = 0 );
		~FitsHeader( void );

   // non-zero codes returned by classify():
   // E end, B blank, * general, H history, C comment, / hanging comment

		operator char*() const;
   char*	operator[] ( int index ) const;
   char		classify( int index = -1 ) const;
   int		length( void ) const;
   int		size( void ) const;

   FitsHeader&	operator= ( const FitsHeader &hdr );
   int		print( int index = -1, FILE *fi = stdout );
   int		extend( int cards );
   int		write( int fd );
   int		read( int fd );

   // codes for seek: ^ first, $ last, > next, < previous
   FitsHeader&	search( const char *str, int seek = '^' );
   FitsHeader&	key( const char *str, int seek = '^' );

   char*	get( void );
   int		operator>> ( char& );
   int		operator>> ( short& );
   int		operator>> ( int& );
   int		operator>> ( long& );
   int		operator>> ( float& );
   int		operator>> ( double& );
   int		operator>> ( char* );

   // conversion and precision are passed through to printf
   FitsHeader&	format( double value, int conversion, int precision );
   FitsHeader&	operator<< ( char );
   FitsHeader&	operator<< ( short );
   FitsHeader&	operator<< ( int );
   FitsHeader&	operator<< ( long );
   FitsHeader&	operator<< ( float );
   FitsHeader&	operator<< ( double );
   FitsHeader&	operator<< ( const char * );

   int	insert( const char *key,const char *rem=0,const char *val=(char*)-1 );
   int	append( const char *key,const char *rem=0,const char *val=(char*)-1 );
   int	update( const char *key,const char *rem=0,const char *val=(char*)-1 );
};

#endif /* _libsite_h */
