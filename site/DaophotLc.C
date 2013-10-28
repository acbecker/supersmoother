/* written by John Doug Reynolds, April 1998 */
/* Normalizer::normalize() written by Mark Pratt, June 1995 */

#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <dirent.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include "libsite.H"

/* --- DaophotLc Extractor ------------------------------------------------- */

class ObjList {
private:
   int			max;
   int			num;
   char			**obj;

public:
   ObjList( void );
   ~ObjList( void );

   int		add( const char *object );
   int		count( void ) const { return num; };
   int		operator[] ( const char *object ) const;
   const char*	operator[] ( int index ) const;
};

class Normalizer {
public:
   int			nobs;		// number of observations
   ObjList		norm;		// list of normalization stars
   int			*toc;		// sorted list of observations
   float		*mags;		// float mags[obs][star]
   float		*errs;		// float errs[obs][star]
   char			*scratch;	// scratch buffer

   inline float&	mag( int obs, int star );   // no bounds checking
   inline float&	err( int obs, int star );   // no bounds checking

   int norm_error;
   int verbose;
   
   //--- Size stuff ---
   int Tobs;		// duplicated from nobs
   int Nref;		// duplicated from norm.count
   int cnobs;		// clipped number of observations
   int cnref;		// clipped number of ref stars

   //--- flags for removing rows, cols and points from mags and err ---
   int *okobs;		// [Tobs]
   int *okstar;		// [Nref]	     
   int **ok;		// [Tobs][Nref] -- nr style matrix (0 offset)

   //--- guess normalization
   float *avgdmag;	// [Tobs] -- mag offset from "best" observation
   float *dmagchi;	// [Tobs] -- chi2 / dof for points defining avgdmag[t]
   float *avgmag;	// [Nref] -- mag of ref stars normalized to best obs
   float *magchi;	// [Nref] -- chi2 / dof for all measurements of star
   float *vbuf;		// [max(Nref,Tobs)] -- sorting vector
   int *nref;		// [Tobs]

   //--- Numerical Recipes interface matrices ie Ax = B ---
   float **A;		// [Tobs][Tobs] -- nr style matrix (1 offset)
   float **B;		// [1][Tobs] -- nr style matrix (1 offset)

   Normalizer( int verbose );
   ~Normalizer( void );

   int		extract( const char *path, DaophotLc &lc );
   int		normalize( DaophotLc &lc );
   int		allocate( void );
   int 		fillAB( void );			
   int 		relaxnorm( const char * );
   int		printok( void );
};

/* --- DaophotLc methods --------------------------------------------------- */

DaophotLc::DaophotLc( void )
{
   nobs = 0;
   object = filter = 0;
   phot = 0;
}

void	DaophotLc::erase( void )
{
   if ( object ) free( object );
   if ( filter ) free( filter );
   if ( phot ) free( (char*) phot );
   errcode = 0;
   nobs = 0;
   object = filter = 0;
   phot = 0;
}

int	DaophotLc::extract( const char *path
			    , const char *obj, const char *fltr, int verbose )
{
   erase();
   if ( !obj  ||  !fltr ) return errcode = errParams;

   object = strdup( obj );
   filter = strdup( fltr );

   Normalizer N( verbose );
   if ( errcode = N.extract( path, *this ) ) return errcode;
   return errcode = N.normalize( *this );
}

int	DaophotLc::print( FILE *fi, int skip, int count, int flags ) const
{
   if ( !fi ) return 1;

   if ( skip >= 0 ) {
      if ( count < 0 ) count = nobs - skip;
      if ( skip + count > nobs ) count = nobs - skip;
   }
   else if ( count < 0 ) count = nobs;

   if ( flags & 1 )
      fprintf( fi, "! %d observations of %s in %s\n", count, object, filter );

   if ( (flags & 3) == 3 ) fprintf( fi, "#\n" );

   if ( flags & 2 )
      fprintf( fi
	       , "#%8s%8s%8s%11s%9s%9s%7s%9s%9s %8s%8s%8s%7s%8s%10s%7s%7s\n"
	       , "date", "mag", "err", "sky", "xpix", "ypix", "amass"
	       , "xFWHM", "yFWHM", "chi2", "sharp", "dMag", "dErr", "dChi2"
	       , "obsid", "norms", "expsr" );

   if ( skip >= 0 )
      for ( int o = skip; o < skip+count; ++o ) {
	 DaophotLc::Photometry &p = phot[o];

	 fprintf( fi,
		  "%9.3f%8.3f%8.3f%11.3f%9.3f%9.3f%7.3f%9.5f%9.5f "
		  "%8.3f%8.3f%8.3f%7.3f%8.3f%10ld%7hd%7hd\n"
		  , machotime(p.date), p.mag, p.err, p.sky, p.xpix, p.ypix
		  , p.airmass, p.xFWHM, p.yFWHM, p.chi2, p.sharpness, p.dMag
		  , p.dMagErr, p.dMagChi2, p.obsid, p.normstars, p.exposure );
      }

   return 0;
}

int	DaophotLc::read( int fd )
{
   erase();
   int nbytes;

   // object
   if ( ::read( fd, (char*) &nbytes, sizeof(int) ) != sizeof(int) )
      return errcode = errFormat;
   if ( !(object = (char*) malloc( nbytes + 1 )) )
      return errcode = errMemory;
   if ( ::read( fd, object, nbytes ) != nbytes )
      return errcode = errFormat;
   object[nbytes] = '\0';

   // filter
   if ( ::read( fd, (char*) &nbytes, sizeof(int) ) != sizeof(int) )
      return errcode = errFormat;
   if ( !(filter = (char*) malloc( nbytes + 1 )) )
      return errcode = errMemory;
   if ( ::read( fd, filter, nbytes ) != nbytes )
      return errcode = errFormat;
   filter[nbytes] = '\0';

   // phot
   if ( ::read( fd, (char*) &nbytes, sizeof(int) ) != sizeof(int) )
      return errcode = errFormat;
   if ( nbytes % sizeof( DaophotLc::Photometry ) )
      return errcode = errFormat;
   if ( !(phot = (DaophotLc::Photometry*) malloc( nbytes )) )
      return errcode = errMemory;
   if ( ::read( fd, (char*) phot, nbytes ) != nbytes )
      return errcode = errFormat;

   nobs = nbytes / sizeof( DaophotLc::Photometry );
   return 0;
}

int	DaophotLc::read( const char *lcdir, const char *obj, const char *fltr )
{
   erase();

   char fname[ MAXPATHLEN+1 ];
   if ( ! fltr ) return errcode = errParams;
   if ( ! lcdir_star_dir(lcdir,obj,'r',fname) ) return errcode = errParams;
   if ( strlen(fname)+strlen(fltr)+7 > MAXPATHLEN ) return errcode = errParams;
   strcat( strcat( fname,  "/daolc." ), fltr );

   int fd = open( fname, O_RDONLY );
   if ( fd == -1 ) return errcode = -errno;
   read( fd );
   close( fd );

   return errcode;
}

int	DaophotLc::write( int fd )
{
   if ( errcode ) return errcode;
   if ( !object || !filter ) return errcode = errFormat;

   int objlen = strlen( object );
   int fltrlen = strlen( filter );
   int nbytes = nobs * sizeof(DaophotLc::Photometry);

   if (
          ::write( fd, (char*) &objlen, sizeof(int) ) != sizeof(int)
       || ::write( fd, object, objlen ) != objlen
       || ::write( fd, (char*) &fltrlen, sizeof(int) ) != sizeof(int)
       || ::write( fd, filter, fltrlen ) != fltrlen
       || ::write( fd, (char*) &nbytes, sizeof(int) ) != sizeof(int)
       || ::write( fd, (char*) phot, nbytes ) != nbytes
      )
      return errcode = -errno;

   return 0;
}

int	DaophotLc::write( const char *lcdir )
{
   if ( errcode ) return errcode;
   if ( !lcdir ) return errcode = errParams;
   if ( !object || !filter || !phot ) return errcode = errFormat;

   char fname[ MAXPATHLEN+1 ];
   if ( !lcdir_star_dir(lcdir,object,'w',fname) ) return errcode = errAccess;
   if ( strlen(fname)+strlen(filter)+7>MAXPATHLEN ) return errcode = errParams;
   strcat( strcat( fname,  "/daolc." ), filter );

   int fd = open( fname, O_RDWR | O_CREAT, 0666 );
   if ( fd == -1 ) return errcode = -errno;
   write( fd );
   close( fd );

   return errcode;
}

/* --- DaophotSet methods -------------------------------------------------- */

DaophotSet::DaophotSet( void )
{
   num = max = 0;
   lc = 0;
}

void	DaophotSet::erase( void )
{
   for ( int i = 0; i < num; ++i ) delete lc[i];
   num = 0;
}

DaophotSet::~DaophotSet( void )
{
   erase();
   delete [] lc;
}

DaophotLc*	DaophotSet::operator[] ( int index ) const
{
   return index >= 0 ? index < num ? lc[index] : 0 : 0;
}

int	DaophotSet::unhook( const DaophotLc *dlc )
{
   for ( int i = 0; i < num; ++i )
      if ( lc[i] == dlc ) {
	 while ( ++i < num ) lc[i-1] = lc[i];
	 lc[--num] = 0;
	 return 0;
      }
   return 1;
}

int	DaophotSet::append( DaophotLc *dlc )
{
   if ( num == max ) {
      DaophotLc **temp = new DaophotLc* [max+10];
      if ( !temp ) return SiteErr::errMemory;
      max += 10;
      if ( lc ) {
	 memcpy( (char*) temp, (char*) lc, num * sizeof( DaophotLc* ) );
	 delete [] lc;
      }
      lc = temp;
   }

   lc[num++] = dlc;
   return 0;
}   

int    DaophotSet::extract( const char *path, const char *object, int verbose )
{
   int error = 0;
   if ( !path || !object ) return SiteErr::errParams;

   DIR *dir;
   struct dirent *d;
   char *scratch = new char [1024];
   sprintf( scratch, "%s/stars/%s", path, object );
   if ( !(dir = opendir(scratch)) ) error = SiteErr::errMissing;

   while ( !error  &&  (d = readdir(dir)) ) {
      struct stat sbuf;
      sprintf( scratch, "%s/stars/%s/%s/normal", path, object, d->d_name );
      if ( stat( scratch, &sbuf ) ) continue;
      DaophotLc *dlc = new DaophotLc;
      if ( !dlc ) error = SiteErr::errMemory;
      if ( !error )
	 if ( dlc->extract(path,object,d->d_name,verbose) || append(dlc) )
	    delete dlc;
   }

   delete [] scratch;
   if ( dir ) closedir( dir );
   return error;
}

int	DaophotSet::read( const char *lcdir, const char *object )
{
   char *path = lcdir_star_dir( lcdir, object );
   if ( !path ) return SiteErr::errParams;

   DIR *dir;
   struct dirent *d;
   if ( !(dir = opendir(path)) ) return SiteErr::errMissing;

   int error = 0;
   while ( !error  &&  (d = readdir(dir)) )
      if ( strncmp( d->d_name, "daolc.", 6 ) == 0 ) {
	 DaophotLc *dlc = new DaophotLc;
	 if ( !dlc ) error = SiteErr::errMemory;
	 if ( !error ) error = dlc->read( lcdir, object, d->d_name+6 );
	 if ( !error ) error = append( dlc );
	 if ( error ) delete dlc;
      }

   closedir( dir );
   delete [] path;
   return error;
}

int	DaophotSet::write( const char *lcdir )
{
   int error = 0;
   for ( int i = 0; !error && i < num; ++i ) error = lc[i]->write( lcdir );
   return error;
}

/* --- ObjList methods ----------------------------------------------------- */

ObjList::ObjList( void )
{
   num = max = 0;
   obj = 0;
}

ObjList::~ObjList( void )
{
   for ( int o = 0; o < num; ++o ) delete [] obj[o];
   delete [] obj;
}

const char*	ObjList::operator[] ( int index ) const
{
   if ( index < 0  ||  index >= num ) return 0;
   return obj[index];
}

int	ObjList::operator[] ( const char *object ) const
{
   for ( int o = 0; o < num; ++o )
      if ( strcmp( object, obj[o] ) == 0 ) return o;
   return -1;
}

int	ObjList::add( const char *object )
{
   int o = (*this)[ object ];
   if ( o != -1 ) return o;
   if ( num == max ) {
      char **temp = new char* [max+100];
      if ( !temp ) return -1;
      max += 100;
      if ( obj ) {
	 memcpy( (char*) temp, (char*) obj, num * sizeof(char*) );
	 delete [] obj;
      }
      obj = temp;
   }
   obj[num] = new char [ strlen(object) + 1 ];
   strcpy( obj[num], object );
   return num++;
}

/* --- Normalizer methods -------------------------------------------------- */

#define SIGCLIP	4

static int	floatcompare( const void *i, const void *j )
{
   float f = *((const float*)i) - *((const float*)j);
   return f >= 0 ? f == 0 ? 0 : 1 : -1;
}
/*--- numerical recipes routines included in source below ---*/
int gaussj(float **a, int n, float **b, int m);
float **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);

Normalizer::Normalizer( int VERBOSE )
{
   verbose = VERBOSE;
   scratch = new char [1024];
   toc = 0;
   mags = errs = 0;
   nobs = 0;
   Tobs = Nref = 0;
   okobs = okstar = nref = 0;
   avgmag = magchi = avgdmag = dmagchi = vbuf = 0;
   ok = 0;
   A = B = 0;
}

int	Normalizer::allocate( void )
{
   int count = nobs * norm.count();
   if ( !count ) return 0;

   mags = new float [ count ];
   errs = new float [ count ];
   toc = new int [ nobs ];

   memset( (char*) toc, 0, nobs * sizeof(int) );
   for ( int i = 0; i < count; ++i ) mags[i] = errs[i] = NaN;

   Tobs = nobs;
   Nref = norm.count();
   okobs = new int [Tobs];
   okstar = new int [Nref];
   nref = new int [Tobs];
   avgmag = new float [Nref];
   magchi = new float [Nref];
   avgdmag = new float [Tobs];
   dmagchi = new float [Tobs];
   vbuf = new float [Nref > Tobs ? Nref : Tobs];

   ok = imatrix( 0, Tobs-1, 0, Nref-1);
   A = matrix(1, Tobs, 1, Tobs);
   B = matrix(1, Tobs, 1, 1);

   return 0;
}

Normalizer::~Normalizer( void )
{
   delete [] toc;
   delete [] mags;
   delete [] errs;
   delete [] scratch;

   delete [] okobs;
   delete [] okstar;
   delete [] nref;
   delete [] avgmag;
   delete [] magchi;
   delete [] avgdmag;
   delete [] dmagchi;
   delete [] vbuf;

   if ( ok ) free_imatrix( ok, 0, Tobs-1, 0, Nref-1);
   if ( A ) free_matrix( A, 1, Tobs, 1, Tobs );
   if ( B ) free_matrix( B, 1, Tobs, 1, 1 );
}

float&	Normalizer::mag( int obs, int star )
{
   return *( mags + obs * norm.count() + star );
}

float&	Normalizer::err( int obs, int star )
{
   return *( errs + obs * norm.count() + star );
}

static int intcompare( const void *i, const void *j )
{
   return *((const int*)i) - *((const int*)j);
}

int	Normalizer::extract( const char *path, DaophotLc &lc )
{
   int obsid;
   if ( !scratch ) return SiteErr::errMemory;

   // build the normalization star table

   sprintf( scratch, "%s/stars/%s/%s/normal", path, lc.object, lc.filter );
   FILE *fi = fopen( scratch, "r" );
   if ( !fi ) return SiteErr::errMissing;
   while ( !feof(fi) ) {
      if ( fscanf( fi, "%s", scratch ) == 1  &&  scratch[0] != '#' )
	 norm.add( scratch );
      fgets(scratch, 1024, fi);
   }
   fclose( fi );

   // count the number of observations present

   DIR *dir;
   struct dirent *d;
   sprintf( scratch, "%s/stars/%s/%s/phot", path, lc.object, lc.filter );
   if ( !(dir = opendir(scratch)) ) return SiteErr::errMissing;

   int nbytes;
   while ( d = readdir(dir) )
      if ( sscanf( d->d_name, "%ld%n", &obsid, &nbytes ) == 1
	   && strlen( d->d_name ) == nbytes )
	 ++nobs;
   rewinddir(dir);

   // allocate space and fill the table of contents

   int error = 0;
   nbytes = nobs * sizeof( DaophotLc::Photometry );
   if ( lc.phot = (DaophotLc::Photometry*) malloc( nbytes ) ) {
      memset( (char*) lc.phot, 0, nbytes );
      error = allocate();
   }
   else error = 1;
   
   while ( !error  &&  lc.nobs < nobs  &&  (d = readdir(dir)) )
      if ( sscanf( d->d_name, "%ld%n", &obsid, &nbytes ) == 1
	   && strlen( d->d_name ) == nbytes )
	 toc[ lc.nobs++ ] = obsid;
   closedir( dir );

   if ( !error )
      qsort( (char*) toc, lc.nobs, sizeof(int), intcompare );

   // read in the photometry

   for ( int o = 0; !error  &&  o < lc.nobs; ++o ) {
      sprintf( scratch, "%s/stars/%s/%s/phot/%ld"
	      , path, lc.object, lc.filter, toc[o] );
      FILE *fi = fopen( scratch, "r" );
      if ( !fi ) return SiteErr::errAccess;

      DaophotLc::Photometry &p = lc.phot[o];
      p.obsid = toc[o];
      fgets( scratch, 1024, fi );
      double date;
      float xHWHM, yHWHM;
      int n = sscanf( scratch, "#%lf%f%f%f%hd"
		      , &date, &p.airmass, &xHWHM, &yHWHM, &p.exposure );
      if ( n < 5 ) error = SiteErr::errFormat;
      p.date = timejulian( date );
      p.xFWHM = 2 * xHWHM;
      p.yFWHM = 2 * yHWHM;

      int found = 0;
      while ( !error  &&  fgets( scratch, 1024, fi ) ) {
	 char obj[80];
	 if ( 1 != sscanf( scratch, "%79s", obj ) ) continue;
	 int s = norm[obj];
	 if ( s != -1 ) {
	    if ( 2 > sscanf( scratch, "%*s%f%f", &mag(o,s), &err(o,s) ) )
	       error = SiteErr::errFormat;
	 }
	 else if ( strcmp( obj, lc.object ) == 0 ) {
	    n = sscanf( scratch, "%*s%f%f%f%f%f%f%f"
			, &p.mag, &p.err, &p.sky, &p.xpix
			, &p.ypix, &p.chi2, &p.sharpness );
	    if ( n < 7 ) error = SiteErr::errFormat;
	    found = 1;
	 }
      }
      if ( !error  &&  !found ) error = SiteErr::errFormat;

      fclose( fi );
   }

   return error;
}

int Normalizer::normalize( DaophotLc &lc )
{
   int m, n, t, dof;
   int T;
   int t0;
   float e, f, g, h;

   norm_error = SiteErr::errNone;
   
   if ( !Tobs  ||  !Nref ) return norm_error;
   
   // Mark valid measurements 
   for( t = 0; t < Tobs; ++t )
      for( n = 0; n < Nref; ++n )
	 ok[t][n] = ( !isnan( err( t, n ) ) && !isnan( mag( t, n ) ) );
   printok();

   // Find best observation
   f = Infinity;
   t0 = -99;
   for( t = 0; t < Tobs; ++t ) {
      g = 0.0;
      for( m = n = 0; n < Nref; ++n )
	 if( ok[t][n] ) {
	    g += err( t, n );
	    ++m;
	 }
      g = m ? g / m : NaN;		// average error
      // printf(t ? "%3d: %8.5f (%d)\n" :
      // "\navgerr\n%3d: %8.5f (%d)\n", t, g, m);
      
      if( g < f ) {
	 f = g;				// new minimum average error
	 t0 = t;			// index of best observation
      }
   }

   // printf( "Min avg err = %8.5f, obs # %d\n", f, t0);
   
   if( f == Infinity )			// No good observations
      return norm_error;

   //----------------------------------------------------------------------
   // Guess at mag normalization and clip outlying measurments.  This is
   // done by first finding all image offsets possible from the "best"
   // observation, constructing a synthesized reference star list and
   // then clipping based on reported errors
   // Start with known mags for reference observation and continue to
   // relaxnorm() until no new observations are added 
   //----------------------------------------------------------------------

   for( n = 0; n < Nref; ++n )
      avgmag[n] = mag(t0,n);
   //printf("Initial reference mags:\n");
   // for( n = 0; n < Nref; ++n )
   //    printf("mag(%2d,%2d): %6.3f\n", t0, n, avgmag[n]);

   while( relaxnorm("median") )			// negative should be error!
      ;
   relaxnorm("median");				// once more to set mags
   if( cnobs < 2 || cnref == 0 )	// useless lightcurve
	 return norm_error;

   // Clip SIGCLIP outliers based on estimated normalization.
   // Remove measurments from futher consideration by resetting ok[n][t]
   // Also, remove measurments from stars and observations not normalized
   // in the last try

   // printok();
   // printf("\nBad obs/star and sigma clipping:\n");
   for( t = 0; t < Tobs; ++t ) 
      for( m = n = 0; n < Nref; ++n ) {

	 // the option exists here to clip on dmagchi[t] and magchi[n]
	 
	 ok[t][n] &= okobs[t] && okstar[n];
	 
	 if( ok[t][n] )			// Sigma clip single points
	    ok[t][n] &= fabs(mag(t,n)-avgdmag[t]-avgmag[n])/err(t,n) < SIGCLIP;
      }

   // printok();
   
   // Some stars and/or observations may have been clipped entirely so
   // relaxnorm until number of normalized obs converges again
   cnobs = cnref = 0;
   while( relaxnorm("average") )	// negative should be error!
      if( cnobs < 2 || cnref == 0 )	// useless lightcurve
	 return norm_error;
   
   // Remove bad rows and columns from ok[t][n]
   // Bad rows should be removed for fillAB, bad cols may not matter
   for( t = 0; t < Tobs; ++t )
      for( n = 0; n < Nref; ++n )
	 ok[t][n] &= ( okobs[t] && okstar[n] );

   // Load nr matrix equation: Ax = B
   if( fillAB() ) return norm_error;

   // Add constraint that one of the dmag[t] = 0
   // This is fishy but the existing system is underconstrained
   if( okobs[t0] ) {			// should always be true
      for( T = t = 0; t < Tobs; ++t )
	 if( okobs[t] )	 {
	    ++T;				// nr indexing
	    if( t == t0 ) A[T][T] += 1; 	// A x0 = 0
	 }
   } else				// worry!
      A[1][1] += 1;			// set first okobs dmag = 0

   // Solve Ax = B
   if( gaussj(A, cnobs, B, 1) )
      return norm_error;

   // Synthesize new reference star mags from fit dmags
   for( T = t = 0; t < Tobs; ++t )
	 avgdmag[t] = okobs[t] ? B[++T][1] : NaN;
   for( n = 0; n < Nref; ++n ) { 
      f = avgmag[n] = 0.0;
      for( t = 0; t < Tobs; ++t ) if( ok[t][n] ) {
	 avgmag[n] += (mag(t,n) - avgdmag[t]) / err(t,n);	// m0 = m - dm
	 f += 1 / err(t,n); 
      }
      avgmag[n] = f ? avgmag[n] / f : NaN;		// error wt'd avg
   }

   // Reference star mags, error in mean and chi2 / dof
   // all calculated from reported photometry errors
   if( verbose ) printf("\nFit: Mags\n");
   for( n = 0; n < Nref; ++n ) {
      g = h = 0.0;
      for( m = t = 0; t < Tobs; ++t ) if( ok[t][n] ) {
	 f = mag(t,n) - avgdmag[t] - avgmag[n];
	 g += 1 / ( err(t,n)*err(t,n) );
	 h += f * f / ( err(t,n)*err(t,n) );
	 ++m;					// number of observations
      }
      vbuf[n] = g > 0.0 ? sqrt( 1 / g ) : NaN;	   // error in mean
      magchi[n] = m > 1 ? sqrt( h / (m-1) ) : NaN; // chi2/dof
      if( verbose && m > 1 ) 
	 printf("avgmag[%2d] = %6.3f  esterr: %6.3f chi2: %6.3f  n: %3d\n",
		n, avgmag[n], vbuf[n], magchi[n], m );
   }

   // Dmag error in mean and chi2/dof
   dof = 0;
   e = 0.0;
   for( t = 0; t < Tobs; ++t ) {
      g = h = 0.0;
      for( m = n = 0; n < Nref; ++n )
	 if( ok[t][n] ) {
	    f = mag(t,n) - avgdmag[t] - avgmag[n];	// measured error
	    g += 1 / ( err(t,n)*err(t,n) );
	    h += f * f / ( err(t,n)*err(t,n) );		// chi2
	    e += f * f / ( err(t,n)*err(t,n) );
	    ++dof;
	    ++m;
	 }
      vbuf[t] = g > 0.0 ? sqrt( 1 / g ) : NaN;	   // estimated dmagerr
      dmagchi[t] = m > 1 ? sqrt( h / (m-1) ) : NaN;// chi2/dof
      
   }

   if( verbose ) {
      printf("\nFit: dMags\n");
      for( t = 0; t < Tobs; ++t ) if( nref[t] )
	 printf("avgdmag[%2d] = %6.3f  esterr: %6.3f  chi2: %6.3f  n: %3d\n",
		t, avgdmag[t], vbuf[t], dmagchi[t], nref[t] );
   }

   dof -= cnobs + cnref;
   if( verbose )
      printf("\nFit chi2/dof = %6.3f  (%.3f,%d)\n", sqrt(e/dof), e, dof);
	 
   // Load dmag vector
   for(t = 0; t < Tobs; ++t )  {
      lc.phot[t].dMag = avgdmag[t];
      lc.phot[t].normstars = nref[t];
      lc.phot[t].dMagErr = vbuf[t];
      lc.phot[t].dMagChi2 = dmagchi[t];
      lc.phot[t].mag -= avgdmag[t];		// overwrite original mags
   }

   return norm_error;
}

int Normalizer::relaxnorm( const char *avgmode ) {
   int lastnobs = cnobs;
   int m, n, t, dof;
   float e, f, g, h;
   int avg_flag = strcmp("average", avgmode) == 0;
   
   cnobs = 0;
   // cnobs is the number of observations that contain at least one
   // normalized reference star.  This is the sum of okobs[t]

   cnref = 0;
   // cnref is the number of reference stars that appear in more
   // than one normalized observation.

   // Normalize all obs which contain normalized reference stars
   for( t = 0; t < Tobs; ++t ) {
      avgdmag[t] = 0.0;				// reset mag offset for obs
      okobs[t] = 0;				// guilty until proven innocent
      f = g = 0.0;
      for( m = n = 0; n < Nref; ++n) 
	 if( ok[t][n] && !isnan(avgmag[n]) ) {
	    vbuf[m] = mag(t,n) - avgmag[n];	// dm = m - m0
	    f += vbuf[m] / err(t,n);
	    g += 1 / err(t,n);
	    // printf("Dmag(%2d,%2d): %6.3f (%5.3f)\n",
	    // t, n, vbuf[m], err(t,n));
	    ++m;
	 }
      if( !m )
	 avgdmag[t] = NaN;			// cannot guess this offset
      else if( avg_flag || m < 3 )		// error wt'd average
	 avgdmag[t] = f / g;
      else {					// take median for >= three
	 qsort( (char *) vbuf, m, sizeof(float), floatcompare );
	 avgdmag[t] = vbuf[m/2];
      }
      if(okobs[t] = (m > 0) ) ++cnobs;		// one hit for okobs
   }

   // Normalize all reference stars in expanded set of normalized obs
   for( n = 0; n < Nref; ++n ) {
      avgmag[n] = 0.0;				// reset reference star mags
      okstar[n] = 0;				// guilty until proven innocent
      f = g = 0.0;
      for( m = t = 0; t < Tobs; ++t )
	 if( ok[t][n] && !isnan(avgdmag[t]) ) {
	    vbuf[m] = mag(t,n) - avgdmag[t];	// m0 = m - dm
	    f += vbuf[m] / err(t,n);
	    g += 1 / err(t,n);
	    ++m;
	 }
      if( !m )
	 avgmag[n] = NaN;			// cannot guess this mag
      else if( avg_flag || m < 3 ) 		// error wt'd average
	 avgmag[n] = f / g;
      else {					// take median for >= three
	 qsort( (char *) vbuf, m, sizeof(float), floatcompare );
	 avgmag[n] = vbuf[m/2];
      }
      okstar[n] = m;				// keep track of hits
   }

   // Calculated magchi[n]
   if( verbose ) printf( "\nrelaxnorm: Mags\n");
   for( n = 0; n < Nref; ++n ) if( !isnan(avgmag[n]) ) {
      f = g = h = 0.0;
      for( m = t = 0; t < Tobs; ++t )		// estimate variability of star
	 if( ok[t][n] && okobs[t] ) {
	    f = mag(t,n) - avgdmag[t] - avgmag[n];
	    h += f * f;
	    f /= err(t,n);
	    g += f * f;
	    ++m;
	 }
      h = m ? sqrt(h / m) : NaN;			// rms error
      magchi[n] = (m > 1) ? sqrt(g / (m - 1)) : NaN;	// chi2 / dof
      if( verbose && m)
	 printf("avgmag[%2d] = %6.3f  rmserr: %6.3f chi2: %6.3f  n: %3d\n",
		n, avgmag[n], h, magchi[n], m );
   }

   // If star appears more than once, it is ok and counts toward cnref
   for( cnref = n = 0; n < Nref; ++n ) 
      if( okstar[n] = (okstar[n] > 1) ) ++cnref;

   // Number of usable reference stars in obs[t]
   for( t = 0; t < Tobs; ++t ) {	
      nref[t] = 0;
      for( n = 0; n < Nref; ++n )
	 if( ok[t][n] && okstar[n] )
	    ++nref[t];
   }
 
   // Calculate dmagchi[t]
   dof = 0;
   e = 0.0;
   if( verbose ) printf( "\nrelaxnorm: dMags\n");
   for( t = 0; t < Tobs; ++t ) if( !isnan(avgdmag[t]) ) {
      f = g = h = 0.0;
      for( n = 0; n < Nref; ++n )		// estimate scatter in dmag
	 if( ok[t][n] && !isnan(avgmag[n]) ) {	// include hit once stars
	    f = mag(t,n) - avgdmag[t] - avgmag[n];
	    h += f * f;
	    f /= err(t,n);
	    g += f * f;
	    e += f * f;
	    ++dof;
	    ++m;
	 }
      h = m ? sqrt(h / m) : NaN;			// rms error
      dmagchi[t] = (m > 1) ? sqrt(g / (m - 1)) : NaN;	// chi2 / dof
      if( verbose )
	 printf("avgdmag[%2d] = %6.3f  rmserr: %6.3f  chi2: %6.3f  n: %3d\n",
		t, avgdmag[t], h, dmagchi[t], nref[t] );
   }

   dof -= cnobs + cnref;
   if( verbose )
      printf("\nrelaxnorm chi2/dof = %6.3f  (%.3f,%d)\n", sqrt(e/dof), e, dof);
      
   //printf("\tstatus: %d obs normalized with %d stars\n", cnobs, cnref );
   //printf("Normalizer::relaxnorm returning %d\n", cnobs-lastnobs);
   
   return cnobs - lastnobs;  	// return increase in number of norm obs
}
   
int Normalizer::fillAB( void ) {
   
   float f;
   int t, i, n;					// std indices
   int T, I;					// nr indices

   //--- Sum vectors ---
   float *Sn = new float [Nref];	// sum over t (1/err2) 
   float *St = new float [Tobs];	// sum over n (1/err2)
   float *Smn = new float [Nref];	// sum over t (mag/err2)
   float *Smt = new float [Tobs];	// sum over n (mag/err2)

   memset( (char *) Sn, 0, Nref * sizeof(float) );
   memset( (char *) St, 0, Tobs * sizeof(float) );
   memset( (char *) Smn, 0, Nref * sizeof(float) );
   memset( (char *) Smt, 0, Tobs * sizeof(float) );

   
   /*--- loop through valid measurments and accumulate sums ---*/
   for( t = 0; t < Tobs; ++t ) 
      for( n = 0; n < Nref; ++n ) 
	 if( ok[t][n] ) {
	    Sn[n] += f = 1 / ( err(t,n)*err(t,n) );
	    St[t] += f;
	    Smt[t] += mag(t,n) * f;
	    Smn[n] += mag(t,n) * f;
	 }

   /*--- loop through all valid measurments and determine vector B ---*/
   for( T = t = 0; t < Tobs; ++t ) {
      if( okobs[t] ) ++T;	 		// T is nr indexing
      else continue;
      B[T][1] = 0.0;				// reset element
      f = 0.0;
      for( n = 0; n < Nref; ++n ) if( ok[t][n] )
	 f += Smn[n] / ( Sn[n] * err(t,n)*err(t,n)  );
      B[T][1] += (f - Smt[t]) / St[t];
   }

   /*--- loop through all valid measurments and determine matrix A ---*/
   for( T = t = 0; t < Tobs; ++t ) {
      if( okobs[t] ) ++T; 			// T is nr indexing (time)
      else continue;
      for( I = i = 0; i < Tobs; ++i ) {	
	 if( okobs[i] ) ++I;			// I is nr indexing (time)
	 else continue;	
	 A[T][I] = 0.0;				// reset element
	 f = 0.0;
	 for( n = 0; n < Nref; ++n ) if( ok[t][n] && ok[i][n] ) 
	    f += 1 / ( Sn[n] * err(t,n)*err(t,n) * err(i,n)*err(i,n) );
	 A[T][I] += f / St[t];
	 if( I == T )
	    A[T][I] -= 1;
      }
   }

   /*------------------- print out A and B --------------------------------
   // Diagnostic
   for( T = t = 0; t < Tobs; ++t ) {
      if( okobs[t] ) ++T; 			// T is nr indexing (time)
      else continue;
      printf("%2d:%2d: ", t, T );
      for( I = i = 0; i < Tobs; ++i ) {	
	 if( okobs[i] ) ++I;			// I is nr indexing (time)
	 else continue;
	 printf("%9.4f ", A[T][I] );
      }
      printf("   %9.4f\n", B[T][1] );
   }
   -----------------------------------------------------------------------*/
   delete [] Sn;
   delete [] St;
   delete [] Smn;
   delete [] Smt;

   return norm_error;
}

int Normalizer::printok( void ) {
   int m = 0;
   int n, t;
   if( !verbose ) return 0;
   printf("\n\t\t\tOK Table\nTime\\Stars\n     ");
   for( n = 0; n < Nref; ++n )
      printf(n < Nref-1 ? "%3d" : "%3d\n", n);
   for( t = 0; t < Tobs; ++t ) {
      printf("%3d: ", t);
      for( n = 0; n < Nref; ++n ) {
	 printf(n < Nref-1 ? "%3d" : "%3d\n", ok[t][n] );
	 m += ok[t][n];
      }
   }
   return m;
}

/*--- numerical recipies Gauss-Jordan elimination ---*/

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

int gaussj(float **a, int n, float **b, int m)
{
   int *indxc,*indxr,*ipiv;
   int i,icol,irow,j,k,l,ll;
   float big,dum,pivinv,temp;

   indxc=ivector(1,n);
   indxr=ivector(1,n);
   ipiv=ivector(1,n);
   for (j=1;j<=n;j++) ipiv[j]=0;
   for (i=1;i<=n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
	 if (ipiv[j] != 1)
	    for (k=1;k<=n;k++) {
	       if (ipiv[k] == 0) {
		  if (fabs(a[j][k]) >= big) {
		     big=fabs(a[j][k]);
		     irow=j;
		     icol=k;
		  }
	       } else if (ipiv[k] > 1) {
		  printf("gaussj: Singular Matrix-1\n");
		  return SiteErr::errMath;
	       }
	    }
      ++(ipiv[icol]);
      if (irow != icol) {
	 for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l]);
	 for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l]);
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (a[icol][icol] == 0.0) {
	 printf("gaussj: Singular Matrix-2\n");
	 return SiteErr::errMath;
      }
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=1;l<=n;l++) a[icol][l] *= pivinv;
      for (l=1;l<=m;l++) b[icol][l] *= pivinv;
      for (ll=1;ll<=n;ll++)
	 if (ll != icol) {
	    dum=a[ll][icol];
	    a[ll][icol]=0.0;
	    for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	    for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
	 }
   }
   for (l=n;l>=1;l--) {
      if (indxr[l] != indxc[l])
	 for (k=1;k<=n;k++)
	    SWAP(a[k][indxr[l]],a[k][indxc[l]]);
   }
   free_ivector(ipiv,1,n);
   free_ivector(indxr,1,n);
   free_ivector(indxc,1,n);
   return SiteErr::errNone;
}
/*--- numerical recipies baggage ---*/
#define NR_END 1
#define FREE_ARG char*

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   float **m;

   /* allocate pointers to rows */
   m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
   if (!m) {
      printf("allocation failure 1 in matrix()\n");
      // norm_error = SiteErr::errMemory;
      return NULL;
   }
   m += NR_END;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
   if (!m[nrl]) {
      printf("allocation failure 2 in matrix()\n");
      // norm_error = SiteErr::errMemory;
      return NULL;
   }
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* return pointer to array of pointers to rows */
   return m;
}
void free_matrix(float **m, long nrl, long , long ncl, long )
/* free a float matrix allocated by matrix() */
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) {
	   printf("allocation failure 1 in matrix()\n");
	   // norm_error = SiteErr::errMemory;
	   return NULL;
	}
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) {
	   printf("allocation failure 2 in matrix()\n");
	   // norm_error = SiteErr::errMemory;
	   return NULL;
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_imatrix(int **m, long nrl, long , long ncl, long )
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
   int *v;

   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) {
      printf( "allocation failure in ivector()\n");
      // norm_error = SiteErr::errMemory;
      return NULL;
   }
   return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long )
/* free an int vector allocated with ivector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}
