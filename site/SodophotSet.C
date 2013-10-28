/* written by John Doug Reynolds, March 1998 */

#include <sys/mman.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include "Macho.h"
#include "libsite.H"

#define BINARY_FORMAT_VERSION 2

/*	Structure of Binary Format Versions 1 and 2
 *	----------------------------------------------------------------------
 *	int			binary format version number
 *	Dimensions		allocated dimensions: nseq, nobs, nchk
 *	int[nchk]		chunk numbers for respective ChunkObs arrays
 *	Header
 *	Star[nseq]
 *	Observation[nobs]
 *	ChunkObs[nchk][nobs]
 *	Photometry[nseq][nobs][2]
 */

struct Star_Version1 {
   unsigned short		seqn;		// star sequence number
   unsigned char		wcid[2];	// canonical chunkid [color]
   SodophotSet::Coordinates	pix[2];		// template coords [color]
   int				ra, dec;	// milli-arcseconds
};

struct ChunkObs_Version1 {
   int				tobsid;		// template observation number
   float			matrix[6];	// coordinate transform matrix
   float			seeing;		// sqrt( PSF area )
   float			avesky;		// average sky level
};

/* --- SodophotSet methods ------------------------------------------------- */

SodophotSet::SodophotSet( void )
{
   bfv = 0; map = 0; hdr = 0; star = 0; obs = 0; phot = 0;
   memset( chnk, 0, sizeof(chnk) );
   max.s = max.o = max.c = 0;
}

SodophotSet::~SodophotSet( void )
{
   close();
}

int	SodophotSet::mapsize( void )
{
   if ( bfv == 1 ) {
      return
	 sizeof(int) + sizeof(Dimensions) + max.c * sizeof(int)
	 + sizeof(Header)
	 + max.s * sizeof(Star_Version1)
	 + max.o * sizeof(Observation)
	 + max.c * max.o * sizeof(ChunkObs_Version1)
	 + 2 * max.s * max.o * sizeof(Photometry);
   }
   else if ( bfv == 2 ) {
      return
	 sizeof(int) + sizeof(Dimensions) + max.c * sizeof(int)
	 + sizeof(Header)
	 + max.s * sizeof(Star)
	 + max.o * sizeof(Observation)
	 + max.c * max.o * sizeof(ChunkObs)
	 + 2 * max.s * max.o * sizeof(Photometry);
   }
   else return 0;
}

// Returns the binary format version number.

int	SodophotSet::version( void ) const
{
   return bfv;
}

// Returns the maximum allocated dimensions of the object.

SodophotSet::Dimensions   SodophotSet::maximum( void ) const
{
   return max;
}

// Lets the object be treated as a standard C 3d array.

SodophotSet::Photometry ( * SodophotSet::operator[] ( int sindex ) const ) [2]
{
   return (Photometry (*)[2]) phot + sindex * max.o;
}

// Deallocates any memory resources used by the object.

void	SodophotSet::close( void )
{
   if ( ! map ) {
      for ( int i = 0; i < 128; ++i ) if ( chnk[i] ) delete [] chnk[i];
      if ( phot ) delete [] phot;
      if ( obs )  delete [] obs;
      if ( star ) delete [] star;
      if ( hdr )  delete hdr;
   }
   else {
      if ( bfv == 1 ) {
	 for ( int i = 0; i < 128; ++i ) if ( chnk[i] ) delete [] chnk[i];
	 if ( star ) delete [] star;
      }
      munmap( map, mapsize() );
   }

   bfv = 0; map = 0; hdr = 0; star = 0; obs = 0; phot = 0;
   memset( chnk, 0, sizeof(chnk) );
   max.s = max.o = max.c = 0;
   errcode = 0;
}

// Dynamically allocates the object to the specified size.

int	SodophotSet::allocate( int nseq, int nobs )
{
   close();

   if ( nseq < 0 || nobs < 0 ) return errcode = errParams;

   hdr = new Header;
   star = new Star [ nseq ];
   obs = new Observation [ nobs ];
   phot = new Photometry [ 2 * nseq * nobs ];

   if ( !hdr || !star || !obs || !phot ) return errcode = errMemory;
   memset( &hdr->count, 0, sizeof(Dimensions) );

   max.s = nseq;
   max.o = nobs;
   max.c = 128;

   return 0;
}

// Allocates chnk[c] if it does not already exist.

int	SodophotSet::allocate( int cid )
{
   if ( cid < 0 || cid > 127 ) return errcode = errParams;
   if ( chnk[cid] ) return 0;

   if ( map && bfv != 1 ) {
      if ( hdr->count.c >= max.c ) return errcode = errMemory;
      ((int*) (map + sizeof(int) + sizeof(Dimensions)))[hdr->count.c] = cid;
      chnk[cid] = &((ChunkObs*) &obs[max.o])[hdr->count.c++ * max.o];
   }
   else {
      chnk[cid] = new ChunkObs [max.o];
      if ( ! chnk[cid] ) return errcode = errMemory;
      ++hdr->count.c;
   }

   return 0;
}

// Memory maps an existing object in the specified file.

static void reverse( void *dst, const void *src, int width, int count = 1 )
{
   char *dp = (char*) dst;
   char *sp = ((char*) src) + width - 1;
   for ( int i = 0; i < count; ++i ) {
      for ( int j = 0; j < width; ++j ) *dp++ = *sp--;
      sp += 2*width;
   }
}

int	SodophotSet::open( const char *filename, char mode )
{
   close();

   if ( !filename || !strchr("rw",mode) ) return errcode = errParams;

   int prot = (mode == 'r') ? PROT_READ : PROT_READ | PROT_WRITE;
   int fd = ::open( filename, mode == 'r' ? O_RDONLY : O_RDWR );
   if ( fd == -1 ) return errcode = -errno;

   if ( read( fd, (char*) &bfv, sizeof(int) ) != sizeof(int) ) {
      bfv = 0;
      errcode = errFormat;
   }
   else if ( bfv < 1 || bfv > BINARY_FORMAT_VERSION ) {
      int temp = bfv;
      reverse( &bfv, &temp, sizeof(int) );
      bfv = ( bfv >= 1 && bfv <= BINARY_FORMAT_VERSION ) ? -bfv : 0;
      errcode = errFormat;
   }
   else if ( bfv == 1 && mode == 'w' ) {
      errcode = errAccess;
   }
   else if (   read( fd, (char*) &max, sizeof(max) ) != sizeof(max)
	    || (map=(char*)mmap(0,mapsize(),prot,MAP_SHARED,fd,0)) == (char*)-1
	    )
      errcode = -errno;
   ::close( fd );
   if ( errcode ) return errcode;

   if ( bfv == 1 ) {
      int j;
      int *label = (int*) (map + sizeof(int) + sizeof(Dimensions));
      hdr = (Header*) &label[max.c];
      Star_Version1 *star_v1 = (Star_Version1*) &hdr[1];
      if ( !(star = new Star [ max.s ]) ) return errcode = errMemory;
      for ( j = 0; j < max.s; ++j ) {
	 Star &d = star[j];
	 Star_Version1 &s = star_v1[j];
	 d.seqn = s.seqn;
	 memcpy( d.wcid, s.wcid, 2 * sizeof(unsigned char) );
	 memcpy( d.pix, s.pix, 2 * sizeof(struct Coordinates) );
	 d.mag[Red] = d.mag[Blue] = INVALID_MAG;
	 d.ra = s.ra;
	 d.dec = s.dec;
      }
      obs = (Observation*) &star_v1[max.s];
      for ( int i = 0; i < max.c; ++i )
	 if ( label[i] >= 0 && label[i] < 128 ) {
	    ChunkObs_Version1 *cobs_v1;
	    cobs_v1 = &((ChunkObs_Version1*) &obs[max.o])[i * max.o];
	    ChunkObs *cobs = chnk[label[i]] = new ChunkObs [ max.o ];
	    if ( !cobs ) return errcode = errMemory;
	    for ( j = 0; j < max.o; ++j ) {
	       ChunkObs &d = cobs[j];
	       ChunkObs_Version1 &s = cobs_v1[j];
	       d.tobsid = s.tobsid;
	       d.ampid = d.ccdflips = -1;
	       d.offset[0] = d.offset[1] = 255;
	       memcpy( d.matrix, s.matrix, 6 * sizeof(float) );
	       d.seeing = s.seeing;
	       d.avesky = s.avesky;
	    }
	 }
      phot = (Photometry*) &((ChunkObs_Version1*) &obs[max.o])[max.c * max.o];
   }
   else {
      int *label = (int*) (map + sizeof(int) + sizeof(Dimensions));
      hdr = (Header*) &label[max.c];
      star = (Star*) &hdr[1];
      obs = (Observation*) &star[max.s];
      for ( int i = 0; i < max.c; ++i )
	 if ( label[i] >= 0 && label[i] < 128 )
	    chnk[label[i]] = &((ChunkObs*) &obs[max.o])[i * max.o];
      phot = (Photometry*) &((ChunkObs*) &obs[max.o])[max.c * max.o];
   }

   return 0;
}

// Constructs and memory maps a new object of the specified maximum size.

int	SodophotSet::open( const char *filename, int nseq, int nobs, int nchk )
{
   close();
   
   if ( !filename || nseq <= 0 || nobs <= 0 || nchk < 0 || nchk > 128 )
      return errcode = errParams;

   int fd = ::open( filename, O_RDWR | O_CREAT, 0666 );
   if ( fd == -1 ) return errcode = -errno;

   Header header;
   int label[128];
   bfv = BINARY_FORMAT_VERSION;
   max.s = nseq; max.o = nobs; max.c = nchk;
   for ( int i = 0; i < nchk; label[i++] = -1 );
   memset( &header, 0, sizeof(Header) );

   if (    ::write( fd, (char*) &bfv, sizeof(int) ) != sizeof(int)
	|| ::write( fd, (char*) &max, sizeof(max) ) != sizeof(max)
	|| ::write( fd, (char*) label, nchk*sizeof(int) ) != nchk*sizeof(int)
	|| ::write( fd, (char*) &header, sizeof(Header) ) != sizeof(Header)
	|| ftruncate( fd, mapsize() )
       )
      errcode = -errno;
   ::close( fd );

   return errcode ? errcode : open( filename, 'w' );
}

// Writes the object to the specified file.

int	SodophotSet::write( const char *filename ) const
{
   if ( !filename || !phot ) return errParams;
   if ( !hdr || !star || !obs ) return errCorrupt;

   // create the list of chunk labels

   int i, label[128];
   struct Dimensions n = hdr->count;
   for ( i = n.c = 0; i < 128; ++i ) if ( chnk[i] ) label[n.c++] = i;
   if ( n.c != hdr->count.c ) return errCorrupt;

   // initialize some useful variables

   int lbytes = n.c * sizeof(int);
   int sbytes = n.s * sizeof(Star);
   int obytes = n.o * sizeof(Observation);
   int cbytes = n.o * sizeof(ChunkObs);
   int pbytes = 2 * n.o * sizeof(Photometry);
   int version = BINARY_FORMAT_VERSION;

   // open the file and write the whole object

   int err = 0, fd = ::open( filename, O_RDWR | O_CREAT | O_TRUNC, 0666 );
   if ( fd == -1 ) return -errno;

   if (    ::write( fd, (char*) &version, sizeof(int) ) != sizeof(int)
	|| ::write( fd, (char*) &n, sizeof(Dimensions) ) != sizeof(Dimensions)
	|| ::write( fd, (char*) label, lbytes ) != lbytes
	|| ::write( fd, (char*) hdr, sizeof(Header) ) != sizeof(Header)
	|| ::write( fd, (char*) star, sbytes ) != sbytes
	|| ::write( fd, (char*) obs, obytes ) != obytes
       )
      err = -errno;

   for ( i = 0; !err && i < n.c; ++i )
      if ( ::write( fd, (char*) chnk[label[i]], cbytes ) != cbytes )
	 err = -errno;

   for ( i = 0; !err && i < n.s; ++i )
      if ( ::write( fd, (char*) (*this)[i], pbytes ) != pbytes )
	 err = -errno;

   ::close( fd );
   return err;
}

// Writes a subset of the object to the specified file.

int	SodophotSet::write( const char *filename
			    , int count, const int *sindex ) const
{
   if ( !filename || count <= 0 || !sindex || !phot ) return errParams;
   if ( !hdr || !star || !obs ) return errCorrupt;

   // determine which chunks are needed while verifying starindex list

   char hit[128];
   struct Dimensions n;
   int i, prev = -1;

   n.o = hdr->count.o;
   memset( hit, 0, 128 );
   for ( n.s = 0; n.s < count; ++n.s ) {
      i = sindex[n.s];
      if ( i <= prev || i < 0 || i >= hdr->count.s ) return errParams;
      Star &s = star[i];
      if ( s.wcid[Red] < 128 ) hit[ s.wcid[Red] ] = 1;
      if ( s.wcid[Blue] < 128 ) hit[ s.wcid[Blue] ] = 1;
      prev = i;
   }

   // create the list of chunk labels

   int label[128];
   for ( n.c = i = 0; i < 128; ++i )
      if ( hit[i] ) {
	 if ( ! chnk[i] ) return errCorrupt;
	 label[n.c++] = i;
      }

   // open the file and write selected lightcurves

   int err = 0, fd = ::open( filename, O_RDWR | O_CREAT | O_TRUNC, 0666 );
   if ( fd == -1 ) return -errno;

   int nbytes = n.c * sizeof(int);
   int version = BINARY_FORMAT_VERSION;
   if (    ::write( fd, (char*) &version, sizeof(int) ) != sizeof(int)
	|| ::write( fd, (char*) &n, sizeof(n) ) != sizeof(n)
	|| ::write( fd, (char*) label, nbytes ) != nbytes
       )
      err = -errno;

   if ( ! err ) {
      Header h = *hdr; h.count = n;
      if ( ::write( fd, (char*) &h, sizeof(h) ) != sizeof(h) ) err = -errno;
   }

   nbytes = sizeof(Star);
   for ( i = 0; !err && i < count; ++i )
      if ( ::write( fd, (char*) &star[sindex[i]], nbytes ) != nbytes )
	 err = -errno;

   nbytes = n.o * sizeof(Observation);
   if ( !err && ::write( fd, (char*) obs, nbytes ) != nbytes )
      err = -errno;

   nbytes = n.o * sizeof(ChunkObs);
   for ( i = 0; !err && i < n.c; ++i )
      if ( ::write( fd, (char*) chnk[label[i]], nbytes ) != nbytes )
	 err = -errno;

   nbytes = 2 * n.o * sizeof(Photometry);
   for ( i = 0; !err && i < count; ++i )
      if ( ::write( fd, (char*) (*this)[sindex[i]], nbytes ) != nbytes )
	 err = -errno;

   ::close( fd );
   return err;
}

// For internal use by printlc() and printso().

static void display( FILE *fi, int width, int precision, double value )
{
   isnan(value)
      ? fprintf( fi, "%*s", width, "-99" )
      : fprintf( fi, "%*.*f", width, precision, value );
}

// For internal use by printlc().

static void printkey( FILE *fi )
{
   fprintf( fi,
	    "!%13s%9s%8s%12s%14s%13s%10s"
	    "%7s%7s%5s%4s%4s%4s%4s%4s%4s%7s%7s%6s%7s%7s%2s "
	    "%7s%7s%5s%4s%4s%4s%4s%4s%4s%7s%7s%6s%7s%7s%2s  \n"
	    , "Date", "Obsid", "Pier", "Exposure", "Checklist", "Airmass", ""
	    , "rMag", "rErr", "rDS", "rTF", "rCP", "rX2", "rMP", "rCR"
	    , "rA", "rXpix", "rYpix", "rSky", "rFWHM", "rTobs", "r"
	    , "bMag", "bErr", "bDS", "bTF", "bCP", "bX2", "bMP", "bCR"
	    , "bA", "bXpix", "bYpix", "bSky", "bFWHM", "bTobs", "b"
	    );
}

// Print an ascii representation of a star-observation, with -99
// substituted for NaN for greater compatibility with supermongo.
// Use -1 for rflag and bflag to apply default cuts, according to
// which a color-observation passes if neither mag nor err is null.

int	SodophotSet::printso( FILE *fi, int starindex, int obsindex
			      , char coords, int rflag, int bflag ) const
{
   if ( !fi || !phot || starindex < 0 || obsindex < 0 )
      return errParams;
   if ( !hdr || !star || !obs )
      return errCorrupt;
   if ( starindex >= hdr->count.s || obsindex >= hdr->count.o )
      return errParams;

   Star &s = star[ starindex ];
   if ( s.wcid[Red] < 128 && !chnk[ s.wcid[Red] ] ) return errCorrupt;
   if ( s.wcid[Blue] < 128 && !chnk[ s.wcid[Blue] ] ) return errCorrupt;

   Observation& o = obs[ obsindex ];

   char checklist[9];
   sprintf( checklist, "0x%.5x", o.checklist );

   ChunkObs bogus;
   memset( &bogus, 0, sizeof(bogus) );
   bogus.avesky = bogus.seeing = NaN;
   bogus.ampid = -1;

   fprintf( fi, "%14.4f%9d%8s%12d%14s"
	    , machotime(o.date), o.obsid, o.pierside == East ? "East" : "West"
	    , o.exposure, checklist );
   display( fi, 13, 4, o.airmass );
   fprintf( fi, "%10s", "" );

   for ( int clr = 0; clr < 2; ++clr ) {
      StarObs p;
      int flag = clr == 0 ? rflag : bflag;
      ChunkObs &c = s.wcid[clr] < 128 ? chnk[s.wcid[clr]][obsindex] : bogus;

      unpack( p, starindex, obsindex, clr, coords );
      if ( flag == -1 ) flag = isnan(p.mag) || isnan(p.err) ? 0 : 1;

      display( fi, 7, 3, p.mag );
      display( fi, 7, p.err < 0 ? 0 : 3, p.err );
      fprintf( fi, "%5d%4d%4d%4d%4d%4d%4d"
	       , p.dsky, p.type, p.crwd, p.chi2, p.mpix, p.cosm, p.ampid );
      display( fi, 7, 1, p.pix.x );
      display( fi, 7, 1, p.pix.y );
      display( fi, 6, 0, c.avesky );
      display( fi, 7, 3, 0.9 * c.seeing );
      fprintf( fi, "%7d%2d ", c.tobsid, flag );
   }

   fprintf( fi, " \n" );
   return 0;
}

// Print an ascii representation of a star-observation, with -99
// substituted for NaN for greater compatibility with supermongo.  If
// starindex is negative all following arguments are ignored and just
// the key of column labels is printed.  If obsindex is negative only
// the lightcurve summary header is printed.  If count is negative it
// is replaced with the number of observations not preceding obsindex,
// or the total number of observations if obsindex is negative.

int	SodophotSet::printlc( FILE *fi, int starindex
			      , char coords, int obsindex, int count ) const
{
   if ( !fi ) return errParams;

   if ( starindex < 0 ) {
      printkey( fi );
      return 0;
   }

   if ( !phot ) return errParams;
   if ( !hdr || !star || !obs ) return errCorrupt;
   if ( starindex >= hdr->count.s ) return errParams;

   if ( obsindex >= 0 ) {
      if ( count < 0 ) count = hdr->count.o - obsindex;
      if ( count < 0 || obsindex + count > hdr->count.o ) return errParams;
   }
   else if ( count < 0 ) count = hdr->count.o;

   Star &s = star[ starindex ];
   if ( s.wcid[Red] < 128 && !chnk[ s.wcid[Red] ] ) return errCorrupt;
   if ( s.wcid[Blue] < 128 && !chnk[ s.wcid[Blue] ] ) return errCorrupt;

   char buffer[32], *source_name = buffer, *coords_name;
   switch ( hdr->source ) {
   case 0: source_name = "an invalid source"; break;
   case 'e': source_name = "the efficiencies simulations"; break;
   case 'g': source_name = "the first-round offline databases"; break;
   case 'n': source_name = "the Notre Dame reductions"; break;
   case 'o': source_name = "the offline databases"; break;
   case 'p': source_name = "the production databases"; break;
   default: sprintf( buffer, "the '%c' reductions", hdr->source ); break;
   }
   switch ( coords ) {
   case 'c': coords_name = "chip"; break;
   case 'a': coords_name = "amp"; break;
   case 'f': coords_name = "fits"; break;
   default: coords_name = "none"; break;
   }

   fprintf( fi, "!  Star %-32s   RA%14s   DEC%14s\n", starid(starindex)
	    , macho_ra_string(s.ra), macho_dec_string(s.dec) );
   fprintf( fi, "#  %d observations from %s\n", count, source_name );
   fprintf( fi, "#  Field inception date: %.4f\n", machotime(hdr->incept) );
   fprintf( fi, "#  Pixel coordinate system: %s\n", coords_name );
   fprintf( fi, "#\n#  %8s%18s%24s%18s\n"
	    , "Template", "magnitude", "chip coordinates", "west-chunk" );

   for ( int clr = 0; clr < 2; ++clr ) {
      fprintf( fi, "#  %6s: ", clr ? "Blue" : "Red" );
      display( fi, 17, 3
	       , s.mag[clr] == INVALID_MAG ? NaN : 0.001 * s.mag[clr] );
      display( fi, 16, 2, s.pix[clr].x );
      fprintf( fi, ", " );
      display( fi, 7, 2, s.pix[clr].y );
      fprintf( fi, "%15d\n", s.wcid[clr] );
   }

   if ( obsindex < 0 ) return 0;

   fprintf( fi, "#\n" );
   printkey( fi );

   for ( int omax = obsindex + count; obsindex < omax; ++obsindex )
      printso( fi, starindex, obsindex, coords );

   return 0;
}

// Unpacks photometry, flags, and coordinates into a more useful form.

void	SodophotSet::unpack( SodophotSet::StarObs &q, int sindex
			     , int oindex, int color, char coords ) const
{
   const Photometry &p = (*this)[sindex][oindex][color];

   q.mag = (p.mag == INVALID_MAG) ? NaN : 0.001 * p.mag;
   q.err = (p.err < 0) ? (p.err==INVALID_SHORT) ? NaN : p.err : 0.001 * p.err;

   const short &flag1 = p.flag[0];
   if ( flag1 != INVALID_SHORT ) {
      q.dsky = ((int) flag1 + 32200) / 100 - 322;
      q.type = flag1 - 100 * q.dsky;
   }
   else q.dsky = q.type = 999;
   
   const short &flag2 = p.flag[1];
   if ( flag2 != INVALID_SHORT ) {
      int tmp = flag2 + 32768;
      q.crwd = tmp >> 8;
      q.chi2 = tmp & 255;
   }
   else q.crwd = q.chi2 = 999;
   
   const short &flag3 = p.flag[2];
   if ( flag3 != INVALID_SHORT ) {
      int tmp = flag3 + 32768;
      q.mpix = tmp >> 8;
      q.cosm = tmp & 255;
   }   
   else q.mpix = q.cosm = 999;

   int wcid = star[sindex].wcid[color];
   q.pix.x = q.pix.y = NaN;
   q.ampid = -1;

   if ( wcid < 128 && chnk[wcid] ) {
      const ChunkObs &c = chnk[wcid][oindex];
      q.ampid = c.ampid;
      if ( ! coords ) return;

      float chipC[2], ampC[2];
      const Coordinates &pix = star[sindex].pix[color];
      chipC[0] = c.matrix[0] * pix.x + c.matrix[1] * pix.y + c.matrix[4];
      chipC[1] = c.matrix[2] * pix.x + c.matrix[3] * pix.y + c.matrix[5];

      if ( c.ampid != -1 && c.ccdflips != -1 ) {
	 int x = c.ccdflips & 1;
	 ampC[0] = ( c.ccdflips & 2 ) ? 2048 - chipC[x] : chipC[x] - 1;
	 ampC[1] = ( c.ccdflips & 4 ) ? 2048 - chipC[1-x] : chipC[1-x] - 1;
	 if ( ampC[0] >= 1024 ) {
	    q.ampid = c.ampid & 0xE;
	    ampC[0] -= 1024;
	 }
	 else if ( ampC[0] < 1024 ) q.ampid = c.ampid | 0x1;
      }
      else ampC[0] = ampC[1] = NaN;

      if ( coords == 'c' ) {
	 q.pix.x = chipC[0];
	 q.pix.y = chipC[1];
      }
      else if ( coords == 'a' ) {
	 q.pix.x = ampC[0];
	 q.pix.y = ampC[1];
      }
      else if ( coords == 'f' && q.ampid == c.ampid
		&& c.offset[0] != 255 && c.offset[1] != 255 ) {
	 q.pix.x = ampC[0] + c.offset[0] + 1;
	 q.pix.y = ampC[1] + c.offset[1] + 1;
      }
   }
}

// Packs photometry and flags, but ignores ampid and coordinates.

void	SodophotSet::pack( const SodophotSet::StarObs &q
			   , int sindex, int oindex, int color )
{
   struct Photometry &p = (*this)[sindex][oindex][color];

   if ( isnan(q.mag) )
      p.mag = INVALID_MAG;
   else
      p.mag = short(1000 * q.mag + (q.mag < 0 ? -0.5 : 0.5));

   if ( isnan(q.err) )
      p.err = INVALID_SHORT;
   else
      p.err = q.err < 0 ? short(q.err) : short(1000 * q.err + 0.5);

   if ( q.dsky == 999 || q.type == 999 )
      p.flag[0] = INVALID_SHORT;
   else
      p.flag[0] = 100 * q.dsky + q.type;

   if ( q.crwd == 999 || q.chi2 == 999 )
      p.flag[1] = INVALID_SHORT;
   else
      p.flag[1] = (q.crwd - 128) * 256 + q.chi2;

   if ( q.mpix == 999 || q.cosm == 999 )
      p.flag[2] = INVALID_SHORT;
   else
      p.flag[2] = (q.mpix - 128) * 256 + q.cosm;
}

// Returns the starindex of the specified star.

int	SodophotSet::find( unsigned short seqn ) const
{
   if ( !hdr || !star ) return -1;

   int i0 = 0;
   int i1 = hdr->count.s - 1;

   while ( i1 >= i0 ) {
      int i = (i0 + i1) / 2;
      int d = star[i].seqn - seqn;
      if ( d < 0 ) i0 = i + 1;
      else if ( d > 0 ) i1 = i - 1;
      else return i;
   }

   return -1;
}

// Returns a proper ascii name for the specified star.

char*	SodophotSet::starid( int sindex, char *object ) const
{
   static char buf[22];

   if ( !hdr || !star || sindex < 0 || sindex >= hdr->count.s ) return 0;
   if ( !object ) object = buf;

   if ( hdr->source == 'g' )
      sprintf( object, "%d.%d"
	       , 1000 * hdr->field + star[sindex].wcid[Red]
	       , (hdr->tile << 15) + star[sindex].seqn );
   else
      sprintf( object, "%d.%d.%d", hdr->field, hdr->tile, star[sindex].seqn );

   return object;
}

// Extracts the current lightcurve from an Offline database.

int	SodophotSet::extract( SodLcMetaDb &db )
{
   close();

   if ( !db.star ) return errcode = errParams;
   if ( allocate( 1, db.obscount() ) ) return errcode;

   if ( db.hdr->iwjsver < 8 ) {
      hdr->source = 'g';
      hdr->tile = db.hdr->tileno[0] - db.star->tile;
   }
   else {
      hdr->source = 'o';
      hdr->tile = db.star->tile;
   }
   hdr->field = db.field;
   hdr->incept = timemacho(0);
   hdr->count.s = 1;

   star->seqn = db.star->seq;
   star->wcid[Red] = db.hdr->chunk;
   star->wcid[Blue] = (ChunkToken(West,db.hdr->chunk)+Blue).id();
   star->pix[Red].x = db.star->x[Red];
   star->pix[Red].y = db.star->y[Red];
   star->pix[Blue].x = db.star->x[Blue];
   star->pix[Blue].y = db.star->y[Blue];
   star->mag[Red] = INVALID_MAG;
   star->mag[Blue] = INVALID_MAG;
   star->ra = db.star->ra;
   star->dec = db.star->dec;

   if ( allocate( star->wcid[Red] ) ) return errcode;
   if ( allocate( star->wcid[Blue] ) ) return errcode;

   for ( int nobs = 0; db.obsindex != -1; db.next() ) {
      SodLcDb::Observation &info = db.obs[ db.obsindex ];
      SodLcDb::Photometry &pht = db.star->phot[ db.obsindex ];
      if ( !info.time  ||  (pht.rerr == -99 && pht.berr == -99) ) continue;

      Observation &o = obs[nobs];
      Photometry &r = (*this)[0][nobs][Red];
      Photometry &b = (*this)[0][nobs][Blue];
      ChunkObs &cr = chnk[star->wcid[Red]][nobs];
      ChunkObs &cb = chnk[star->wcid[Blue]][nobs];

      o.date = timemacho( info.time );
      o.obsid = info.obsno;
      o.checklist = FlatFielded;
      o.exposure = INVALID_SHORT;
      o.pierside = info.ipier ? East : West;
      o.airmass = info.airmass;

      if ( pht.rerr != -99 ) {
	 r.mag = pht.rmag;
	 r.err = pht.rerr;
	 r.flag[0] = pht.rflag1;
	 r.flag[1] = pht.rflag2;
	 r.flag[2] = pht.rflag3;
      }
      else {
	 r.mag = INVALID_MAG;
	 r.err = r.flag[0] = r.flag[1] = r.flag[2] = INVALID_SHORT;
      }

      if ( pht.berr != -99 ) {
	 b.mag = pht.bmag;
	 b.err = pht.berr;
	 b.flag[0] = pht.bflag1;
	 b.flag[1] = pht.bflag2;
	 b.flag[2] = pht.bflag3;
      }
      else {
	 b.mag = INVALID_MAG;
	 b.err = b.flag[0] = b.flag[1] = b.flag[2] = INVALID_SHORT;
      }

      cr.tobsid = (int) info.tempobs[Red];
      cb.tobsid = (int) info.tempobs[Blue];
      cr.seeing = info.seeing[Red];
      cb.seeing = info.seeing[Blue];
      cr.avesky = info.avsky[Red];
      cb.avesky = info.avsky[Blue];

      cr.ampid = cb.ampid = cr.ccdflips = cb.ccdflips = -1;
      cr.offset[0] = cr.offset[1] = cb.offset[0] = cb.offset[1] = 255;
      for ( int i = 0; i < 6; ++i ) cr.matrix[i] = cb.matrix[i] = NaN;

      hdr->count.o = ++nobs;
   }

   return errNone;
}

#ifndef __production__

int	SodophotSet::extract( ExtractLcOpts, int, const ushort* )
{
   close();
   return errcode = errVersion;
}

int	SodophotSet::extract( ExtractLcOpts, int, int, int )
{
   close();
   return errcode = errVersion;
}

#endif

/* --- miscellaneous support routines -------------------------------------- */

ExtractLcOpts::ExtractLcOpts( void )
{
   active = cache = vmfile = 0;
   minobs = maxobs = quota = 0;
   merge = unzip = fetch = 0;
   field = tile = 0;
   tolerant = 0;
   verbose = 0;
}

// in this version of merge, s1 is assumed to be more recent than s2
// and both sets must contain exactly the same star, and only one star

int merge( SodophotSet &sum, const SodophotSet &s1, const SodophotSet &s2 )
{
   if (    s1.hdr->field != s2.hdr->field || s1.hdr->tile != s2.hdr->tile
	|| s1.hdr->source != s2.hdr->source || s1.hdr->count.s != 1
	|| s2.hdr->count.s != 1 || s1.star->seqn != s2.star->seqn
      )
      return SiteErr::errCorrupt;

   // allocate space in sum to just fit s1 + s2

   int nobs, o1=0, o2=0;
   int nobs1 = s1.hdr->count.o;
   int nobs2 = s2.hdr->count.o;
   int rwcid = s1.star->wcid[Red];
   int bwcid = s1.star->wcid[Blue];

   for ( nobs = 0; o1 < nobs1 && o2 < nobs2; ++nobs ) {
      int d = s1.obs[o1].obsid - s2.obs[o2].obsid;
      if ( d == 0 ) { ++o1; ++o2; }
      else if ( d < 0 ) ++o1;
      else ++o2;
   }
   if ( o1 < nobs1 ) nobs += nobs1 - o1;
   if ( o2 < nobs2 ) nobs += nobs2 - o2;

   if ( sum.allocate( 1, nobs ) ) return sum.errcode;
   if ( rwcid < 128 && sum.allocate( rwcid ) ) return sum.errcode;
   if ( bwcid < 128 && sum.allocate( bwcid ) ) return sum.errcode;

   // merge s1 and s2 into sum

   sum.hdr->field = s1.hdr->field;
   sum.hdr->tile = s1.hdr->tile;
   sum.hdr->count = s1.hdr->count;
   sum.hdr->source = s1.hdr->source;
   sum.hdr->incept = s1.hdr->incept;
   sum.hdr->count.o = nobs;
   sum.star[0] = s1.star[0];

   SodophotSet::ChunkObs bogus;
   bogus.tobsid = 0;
   bogus.ampid = bogus.ccdflips = -1;
   bogus.offset[0] = bogus.offset[1] = 255;
   for ( int i = 0; i < 6; ++i ) bogus.matrix[i] = NaN;
   bogus.seeing = bogus.avesky = NaN;

   for ( int o = o1 = o2 = 0; o < nobs; ++o ) {
      int d, o0;
      if ( o1 < nobs1 && o2 < nobs2 ) {
	 d = s1.obs[o1].obsid - s2.obs[o2].obsid;
	 if ( d == 0 ) { d = -1; ++o2; }
      }
      else d = o1 < nobs1 ? -1 : +1;
      o0 = d < 0 ? o1++ : o2++;
      const SodophotSet &s0 = *( d < 0 ? &s1 : &s2 );

      sum.obs[o] = s0.obs[o0];
      sum[0][o][Red] = s0[0][o0][Red];
      sum[0][o][Blue] = s0[0][o0][Blue];

      int rwcid0 = s0.star->wcid[Red];
      int bwcid0 = s0.star->wcid[Blue];

      if ( rwcid < 128 )
	 if ( rwcid0 < 128 ) {
	    SodophotSet::ChunkObs &c = sum.chnk[rwcid][o];
	    c = s0.chnk[rwcid0][o0];
	    if ( rwcid != rwcid0 )
	       memcpy( c.matrix, bogus.matrix, sizeof(c.matrix) );
	 }
	 else sum.chnk[rwcid][o] = bogus;

      if ( bwcid < 128 )
	 if ( bwcid0 < 128 ) {
	    SodophotSet::ChunkObs &c = sum.chnk[bwcid][o];
	    c = s0.chnk[bwcid0][o0];
	    if ( bwcid != bwcid0 )
	       memcpy( c.matrix, bogus.matrix, sizeof(c.matrix) );
	 }
	 else sum.chnk[bwcid][o] = bogus;
   }

   return 0;
}

// create a new SodophotSet file as a byte-order-reversed copy of the
// original

int sodset_byte_reverse(const char *dstfile, const char *srcfile, int verbose)
{
   SiteErr err;
   const char *context = "sodset_byte_reverse";

   // open the source file

   if ( !dstfile || !srcfile ) {
      err = SiteErr::errParams;
      return verbose ? err.perror(context) : int(err);
   }

   int src = open( srcfile, O_RDONLY );
   if ( src == -1 ) {
      err = -errno;
      return verbose ? err.perror(srcfile) : int(err);
   }

   // determine the binary format version and byte-order of the source file
   // w and m are set such that data read into [w] and reversed into [m] is
   // always arranged so that [0] is in the byte-order of this architecture
   // o and r are associated with w and m respectively, and refer to
   // obverse and reverse

   int bfv[2], w, m;
   if ( read( src, (char*) bfv, sizeof(int) ) != sizeof(int) ) {
      close( src );
      err = SiteErr::errFormat;
      return verbose ? err.perror(srcfile) : int(err);
   }
   if ( *bfv < 1 || *bfv > BINARY_FORMAT_VERSION ) {
      bfv[1] = *bfv;
      reverse( bfv, bfv+1, sizeof(int) );
      if ( *bfv < 1 || *bfv > BINARY_FORMAT_VERSION ) {
	 close( src );
	 err = SiteErr::errFormat;
	 return verbose ? err.perror(srcfile) : int(err);
      }
      w = 1;
   }
   else {
      reverse( bfv+1, bfv, sizeof(int) );
      w = 0;
   }
   m = 1 - w;

   if ( *bfv != 2 ) {
      close( src );
      err = SiteErr::errVersion;
      return verbose ? err.perror(context) : int(err);
   }

   // open the output file and write the binary format version

   int dst = open( dstfile, O_RDWR | O_CREAT, 0666 );
   if ( dst == -1 ) {
      err = -errno;
      close( src );
      return verbose ? err.perror(dstfile) : int(err);
   }

   if ( write( dst, (char*)(bfv+m), sizeof(int) ) != sizeof(int) ) {
      err = -errno;
      if ( verbose ) err.perror(dstfile);
   }

   // read the dimensions, chunk labels, and header

   int label[2][128];
   SodophotSet::Dimensions n;
   SodophotSet::Dimensions max[2];
   SodophotSet::Header hdr[2];

   if ( ! err ) {
      int nbytes = sizeof(SodophotSet::Dimensions);
      if ( read( src, (char*)(max+w), nbytes ) == nbytes ) {
	 reverse( &max[m].s, &max[w].s, sizeof(int) );
	 reverse( &max[m].o, &max[w].o, sizeof(int) );
	 reverse( &max[m].c, &max[w].c, sizeof(int) );
      }
      else err = -errno;
   }

   if ( ! err ) {
      int nbytes = max->c * sizeof(int);
      if ( read( src, (char*)(label+w), nbytes ) == nbytes )
	 reverse( label[m], label[w], sizeof(int), max->c );
      else err = -errno;
   }

   if ( ! err ) {
      int nbytes = sizeof(SodophotSet::Header);
      if ( read( src, (char*)(hdr+w), nbytes ) == nbytes ) {
	 SodophotSet::Header &o = hdr[w];
	 SodophotSet::Header &r = hdr[m];
	 reverse( &r.field, &o.field, sizeof(short) );
	 reverse( &r.tile, &o.tile, sizeof(unsigned short) );
	 reverse( &r.count.s, &o.count.s, sizeof(int) );
	 reverse( &r.count.o, &o.count.o, sizeof(int) );
	 reverse( &r.count.c, &o.count.c, sizeof(int) );
	 reverse( &r.source, &o.source, sizeof(int) );
	 reverse( &r.incept, &o.incept, sizeof(int) );
	 n = hdr->count;
      }
      else err = -errno;
   }

   if ( err && verbose ) err.perror(srcfile);

   // write the dimensions (minimized), chunk labels, and header

   if ( ! err ) {
      int lbytes = n.c * sizeof(int);
      int dbytes = sizeof(SodophotSet::Dimensions);
      int hbytes = sizeof(SodophotSet::Header);

      if ( write( dst, (char*)&hdr[m].count, dbytes ) != dbytes ||
	   write( dst, (char*)(label+m), lbytes ) != lbytes ||
	   write( dst, (char*)(hdr+m), hbytes ) != hbytes )
	 {
	    err = -errno;
	    if ( verbose ) err.perror(dstfile);
	 }
   }

   // read and write the star records

   if ( ! err ) {
      SodophotSet::Star o, r;
      int nbytes = sizeof(SodophotSet::Star);
      for ( int i = 0; !err && i < n.s; ++i )
	 if ( read( src, (char*)&o, nbytes ) == nbytes ) {
	    reverse( &r.seqn, &o.seqn, sizeof(unsigned short) );
	    memcpy( r.wcid, o.wcid, 2 );
	    reverse( &r.pix[0].x, &o.pix[0].x, sizeof(float) );
	    reverse( &r.pix[0].y, &o.pix[0].y, sizeof(float) );
	    reverse( &r.pix[1].x, &o.pix[1].x, sizeof(float) );
	    reverse( &r.pix[1].y, &o.pix[1].y, sizeof(float) );
	    reverse( r.mag, o.mag, sizeof(short), 2 );
	    reverse( &r.ra, &o.ra, sizeof(int) );
	    reverse( &r.dec, &o.dec, sizeof(int) );
	    if ( write( dst, (char*)&r, nbytes ) != nbytes ) {
	       err = -errno;
	       if ( verbose ) err.perror(dstfile);
	    }
	 }
	 else {
	    err = -errno;
	    if ( verbose ) err.perror(srcfile);
	 }
      int skip = nbytes * (max->s - n.s);
      if ( !err && skip > 0 && lseek(src,skip,SEEK_CUR) == -1 ) {
	 err = -errno;
	 if ( verbose ) err.perror(srcfile);
      }
   }

   // read and write the observation records

   if ( ! err ) {
      SodophotSet::Observation o, r;
      int nbytes = sizeof(SodophotSet::Observation);
      for ( int i = 0; !err && i < n.o; ++i )
	 if ( read( src, (char*)&o, nbytes ) == nbytes ) {
	    reverse( &r.date, &o.date, sizeof(int) );
	    reverse( &r.obsid, &o.obsid, sizeof(int) );
	    reverse( &r.checklist, &o.checklist, sizeof(int) );
	    reverse( &r.exposure, &o.exposure, sizeof(short) );
	    reverse( &r.pierside, &o.pierside, sizeof(short) );
	    reverse( &r.airmass, &o.airmass, sizeof(float) );
	    if ( write( dst, (char*)&r, nbytes ) != nbytes ) {
	       err = -errno;
	       if ( verbose ) err.perror(dstfile);
	    }
	 }
	 else {
	    err = -errno;
	    if ( verbose ) err.perror(srcfile);
	 }
      int skip = nbytes * (max->o - n.o);
      if ( !err && skip > 0 && lseek(src,skip,SEEK_CUR) == -1 ) {
	 err = -errno;
	 if ( verbose ) err.perror(srcfile);
      }
   }

   // read and write the chunk records

   if ( ! err ) {
      SodophotSet::ChunkObs o, r;
      int nbytes = sizeof(SodophotSet::ChunkObs);
      int skip = nbytes * (max->o - n.o);
      for ( int i = 0; !err && i < n.c; ++i ) {
	 for ( int j = 0; !err && j < n.o; ++j )
	    if ( read( src, (char*)&o, nbytes ) == nbytes ) {
	       reverse( &r.tobsid, &o.tobsid, sizeof(int) );
	       r.ampid = o.ampid;
	       r.ccdflips = o.ccdflips;
	       memcpy( r.offset, o.offset, 2 );
	       reverse( r.matrix, o.matrix, sizeof(float), 6 );
	       reverse( &r.seeing, &o.seeing, sizeof(float) );
	       reverse( &r.avesky, &o.avesky, sizeof(float) );
	       if ( write( dst, (char*)&r, nbytes ) != nbytes ) {
		  err = -errno;
		  if ( verbose ) err.perror(dstfile);
	       }
	    }
	    else {
	       err = -errno;
	       if ( verbose ) err.perror(srcfile);
	    }
	 if ( !err && skip > 0 && lseek(src,skip,SEEK_CUR) == -1 ) {
	    err = -errno;
	    if ( verbose ) err.perror(srcfile);
	 }
      }
      skip = nbytes * max->o * (max->c - n.c);
      if ( !err && skip > 0 && lseek(src,skip,SEEK_CUR) == -1 ) {
	 err = -errno;
	 if ( verbose ) err.perror(srcfile);
      }
   }

   // read and write the photometry records

   if ( ! err ) {
      SodophotSet::Photometry o, r;
      int nbytes = sizeof(SodophotSet::Photometry);
      int skip = nbytes * (max->o - n.o);
      for ( int i = 0; !err && i < n.s; ++i ) {
	 for ( int j = 0; !err && j < n.o; ++j )
	    for ( int k = 0; !err && k < 2; ++k )
	       if ( read( src, (char*)&o, nbytes ) == nbytes ) {
		  reverse( &r.mag, &o.mag, sizeof(short) );
		  reverse( &r.err, &o.err, sizeof(short) );
		  reverse( r.flag, o.flag, sizeof(short), 3 );
		  if ( write( dst, (char*)&r, nbytes ) != nbytes ) {
		     err = -errno;
		     if ( verbose ) err.perror(dstfile);
		  }
	       }
	       else {
		  err = -errno;
		  if ( verbose ) err.perror(srcfile);
	       }
	 if ( !err && skip > 0 && lseek(src,skip,SEEK_CUR) == -1 ) {
	    err = -errno;
	    if ( verbose ) err.perror(srcfile);
	 }
      }
      skip = 2 * nbytes * max->o * (max->s - n.s);
      if ( !err && skip > 0 && lseek(src,skip,SEEK_CUR) == -1 ) {
	 err = -errno;
	 if ( verbose ) err.perror(srcfile);
      }
   }

   // done

   close( dst );
   close( src );
   return err;
}

// create a new SodophotSet from selected stars and observations
// of the original

int excerpt( SodophotSet &dst, const SodophotSet &src
	     , int scount, const int *sindex, int ocount, const int *oindex )
{
   dst.close();

   if ( !src[0] || (sindex && scount < 0) || (oindex && ocount < 0) )
      return dst.errcode = SiteErr::errParams;
   if ( !src.hdr || !src.star || !src.obs )
      return dst.errcode = SiteErr::errCorrupt;

   // create lists of stars and observations if necessary

   int *starindex = 0, *obsindex = 0;
   if ( ! sindex ) {
      scount = src.hdr->count.s;
      sindex = starindex = new int [scount];
      if ( sindex ) for ( int i = 0; i < scount; ++i ) starindex[i] = i;
   }
   if ( ! oindex ) {
      ocount = src.hdr->count.o;
      oindex = obsindex = new int [ocount];
      if ( oindex ) for ( int i = 0; i < ocount; ++i ) obsindex[i] = i;
   }
   if ( !sindex || !oindex ) dst.errcode = SiteErr::errMemory;

   // determine which chunks are needed and verify lists

   char hit[128];
   int h, i, prev;
   memset( hit, 0, 128 );

   if ( ! dst.errcode )
      for ( prev = -1, h = 0; h < scount; ++h )
	 if ( (i=sindex[h]) <= prev || i < 0 || i >= src.hdr->count.s ) {
	    dst.errcode = SiteErr::errParams;
	    break;
	 }
	 else {
	    SodophotSet::Star &s = src.star[prev=i];
	    if ( s.wcid[Red] < 128 ) hit[ s.wcid[Red] ] = 1;
	    if ( s.wcid[Blue] < 128 ) hit[ s.wcid[Blue] ] = 1;
	 }

   if ( ! dst.errcode )
      for ( prev = -1, h = 0; h < ocount; ++h )
	 if ( (i=oindex[h]) <= prev || i < 0 || i >= src.hdr->count.o ) {
	    dst.errcode = SiteErr::errParams;
	    break;
	 }
	 else prev = i;

   // allocate dst and verify src

   if ( !dst.errcode && !dst.allocate(scount,ocount) )
      for ( i = 0; i < 128; ++i )
	 if ( hit[i] ) {
	    if ( ! src.chnk[i] ) {
	       dst.errcode = SiteErr::errCorrupt;
	       break;
	    }
	    if ( dst.allocate(i) ) break;
	 }

   // copy from src to dst

   if ( ! dst.errcode ) {
      int ccount = dst.hdr->count.c;
      *dst.hdr = *src.hdr;
      dst.hdr->count.s = scount;
      dst.hdr->count.o = ocount;
      dst.hdr->count.c = ccount;

      for ( i = 0; i < scount; ++i )
	 dst.star[i] = src.star[ sindex[i] ];

      for ( i = 0; i < ocount; ++i )
	 dst.obs[i] = src.obs[ oindex[i] ];

      for ( h = 0; h < 128; ++h )
	 if ( hit[h] )
	    for ( i = 0; i < ocount; ++i )
	       dst.chnk[h][i] = src.chnk[h][ oindex[i] ];

      for ( h = 0; h < scount; ++h )
	 for ( i = 0; i < ocount; ++i )
	    memcpy( dst[h][i], src[sindex[h]][oindex[i]]
		    , 2*sizeof(SodophotSet::Photometry) );
   }

   // cleanup

   if ( starindex ) delete [] starindex;
   if ( obsindex ) delete [] obsindex;
   return dst.errcode;
}
