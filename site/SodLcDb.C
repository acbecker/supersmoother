/* written by John Doug Reynolds, April 1996 */

#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include "libsite.H"

/* --- SodLcIndex methods -------------------------------------------------- */

SodLcIndex::SodLcIndex( const char *path, const char *name )
{
   map = 0;
   open( path, name );
}

int	SodLcIndex::close( void )
{
   if ( map ) munmap( map, size );
   map = 0;
   return errcode = 0;
}

int	SodLcIndex::open( const char *path, const char *name )
{
   close();
   if ( !path ) return errcode = errParams;
   if ( !name ) name = "Lightcurve/index";
   char fname[1024];
   sprintf( fname, "%s/%s", path, name );
   struct stat statbuf;
   if ( stat( fname, &statbuf ) ) return errcode = -errno;
   size = (int) statbuf.st_size;
   int fd = ::open( fname, O_RDONLY );
   if ( fd == -1 ) return errcode = -errno;
   map = (char*) mmap( NULL, size, PROT_READ, MAP_SHARED, fd, 0 );
   ::close( fd );
   if ( map == (char*) -1 ) return errcode = -errno;
}

int	SodLcIndex::count( int tile ) const
{
   if ( tile < 0  ||  tile >= MaxTile ) return 0;
   return ((int*)map)[tile] & 0xFF;
}   

SodLcIndex::LcIndex*	SodLcIndex::operator[] ( int tile ) const
{
   if ( tile < 0  ||  tile >= MaxTile ) return 0;
   int address = ((int*)map)[tile];
   LcIndex *index0 = (LcIndex*) ((int*)map + MaxTile);
   return address ? index0 + (address >> 8) : 0;
}

/* --- SodLcDb methods ------------------------------------------------ */

int	SodLcDb::close( void )
{
   if ( hdr ) munmap( (char*) hdr, filesize );
   delete [] fname;
   fname = 0;
   hdr = 0;
   obs = 0;
   lc0 = 0;
   sizeof_Lightcurve = 0;
   return errcode = 0;
}

SodLcDb::SodLcDb( void )
{
   hdr = 0;
   fname = 0;
   close();
}

SodLcDb::SodLcDb( const char *FNAME )
{
   hdr = 0;
   fname = 0;
   open( FNAME );
}

SodLcDb::SodLcDb( const char *FNAME, int nstars, int nobs )
{
   hdr = 0;
   fname = 0;
   open( FNAME, nstars, nobs );
}

int	SodLcDb::open( const char *FNAME )
{
   close();

   fname = new char [ strlen(FNAME) + 1 ];
   if ( !fname ) return errcode = errMemory;
   strcpy( fname, FNAME );

   struct stat statbuf;
   if ( stat( fname, &statbuf ) ) return errcode = -errno;
   filesize = (int) statbuf.st_size;
   int fd = ::open( fname, O_RDONLY );
   if ( fd == -1 ) return errcode = -errno;
   hdr = (Header*) mmap( NULL, filesize, PROT_READ, MAP_SHARED, fd, 0 );
   ::close( fd );
   if ( hdr == (Header*) -1 ) return errcode = -errno;

   sizeof_Lightcurve =
      sizeof( Lightcurve ) + int(hdr->nobs-1) * sizeof( SodLcDb::Photometry );

   obs = (Observation*) &hdr[1];
   lc0 = (char*) &obs[hdr->nobs];

   return errcode;
}

int	SodLcDb::open( const char *FNAME, int nstars, int nobs )
{
   close();

   if ( nstars < 1  ||  nobs < 1 ) return errcode = errParams;

   fname = new char [ strlen(FNAME) + 1 ];
   if ( !fname ) return errcode = errMemory;
   strcpy( fname, FNAME );

   sizeof_Lightcurve =
      sizeof( SodLcDb::Lightcurve ) + (nobs-1) * sizeof( SodLcDb::Photometry );

   filesize =
      sizeof( SodLcDb::Header )
      + nobs * sizeof( SodLcDb::Observation )
      + nstars * sizeof_Lightcurve;

   int fd = ::open( fname, O_RDWR | O_CREAT, 0666 );
   if ( fd == -1  ||  ftruncate( fd, filesize ) ) return errcode = -errno;
   int protections = PROT_READ | PROT_WRITE;
   hdr = (Header*) mmap( NULL, filesize, protections, MAP_SHARED, fd, 0 );
   ::close( fd );
   if ( hdr == (Header*) -1 ) return errcode = -errno;

   obs = (Observation*) &hdr[1];
   lc0 = (char*) &obs[nobs];
   hdr->nstars = nstars;
   hdr->nobs = nobs;

   hdr->_rec1_start = hdr->_rec1_end
      = (char*) &hdr->_rec1_end - (char*) &hdr->_rec1_start - 4;
   hdr->_rec2_start = hdr->_rec2_end
      = (char*) &hdr->_rec2_end - (char*) &hdr->_rec2_start - 4;
   hdr->_rec3_start = hdr->_rec3_end
      = (char*) &hdr->_rec3_end - (char*) &hdr->_rec3_start - 4;
   hdr->_rec4_start = hdr->_rec4_end
      = (char*) &hdr->_rec4_end - (char*) &hdr->_rec4_start - 4;
   hdr->_rec5_start = hdr->_rec5_end
      = (char*) &hdr->_rec5_end - (char*) &hdr->_rec5_start - 4;
   hdr->_rec6_start = hdr->_rec6_end
      = (char*) &hdr->_rec6_end - (char*) &hdr->_rec6_start - 4;
   hdr->_rec7_start = hdr->_rec7_end
      = (char*) &hdr->_rec7_end - (char*) &hdr->_rec7_start - 4;
   int i;
   for ( i = 0; i < nobs; ++i )
      obs[i]._rec_start = obs[i]._rec_end = sizeof( SodLcDb::Observation ) - 8;
   for ( i = 0; i < nstars; ++i ) {
      SodLcDb::Lightcurve& lc = *operator[](i);
      lc._rec1_start = lc._rec1_end
	 = (char*) &lc._rec1_end - (char*) &lc._rec1_start - 4;
      lc._rec2_start = lc._rec2_end
	 = (char*) &lc._rec2_end - (char*) &lc._rec2_start - 4;
      lc._rec3_start = lc._rec3_end
	 = (char*) &lc._rec3_end - (char*) &lc._rec3_start - 4;
      lc._rec4_start = nobs * sizeof( SodLcDb::Photometry );
      *( (int*) &lc.phot[nobs] ) = lc._rec4_start;
   }

   return errcode;
}

// binary search
SodLcDb::Lightcurve* SodLcDb::find( int tile, int seq, int index, int count )
{
   if ( !hdr ) return 0;

   int min = index;
   int max = index + count - 1;
   int sgn = hdr->iwjsver < 8 ? -1 : 1;

   if ( min < 0 ) min = 0;
   if ( max >= hdr->nstars ) max = hdr->nstars - 1;

   while ( max >= min ) {
      int i = (max + min) / 2;
      Lightcurve *lc = (*this)[i];
      int diff = tile == lc->tile ? lc->seq - seq : sgn * (lc->tile - tile);
      if ( diff == 0 )  return lc;
      else if ( diff < 0 ) min = i + 1;
      else max = i - 1;
   }

   return 0;
}

SodLcDb::Lightcurve*	SodLcDb::find( int tile, int seq )
{
   if ( !hdr ) return 0;
   return find( tile, seq, 0, hdr->nstars );
}

/* --- SodLcMetaDb methods ------------------------------------------------- */

SodLcMetaDb::SodLcMetaDb( void )
{
   db = 0; toc = 0;
   close();
}

SodLcMetaDb::~SodLcMetaDb( void )
{
   close();
}

int	SodLcMetaDb::close( void )
{
   if ( db ) delete [] db;
   if ( toc ) delete [] toc;
   index = nobs = ndb = 0;
   obsindex = primary = segment = -1;
   field = 0; fname = 0; hdr = 0; obs = 0;
   db = 0; toc = 0; star = 0;
   return 0;
}

static int	select_db_file( const char *name, int field, int chunk )
{
   int p, f, c, n;
   const char *suffix = name + strlen(name);
   while ( suffix > name && *--suffix != '_' );
   while ( suffix > name && *--suffix != '_' );
   p = sscanf( suffix, "_f%d_c%d%n", &f, &c, &n );
   return (p == 2 && f == field && c == chunk && suffix[n] == '\0');
}

/*
  SodLcMetaDb::open()
     wjstype: 1 == Round1 (iwjsver < 8), 2 == Round2+, -1 == first type seen
*/

int	SodLcMetaDb::open(const char *path, int fld, int chunk, int wjstype)
{
   close();

   // count the number of database segments and allocate space for them

   DIR *dir;
   struct dirent *d;
   char scratch[512];
   sprintf( scratch, "%s/Lightcurve/f_%d", path, field = fld );
   if ( !(dir = opendir(scratch)) ) return errcode = -errno;
   int maxd=0;
   while ( d = readdir(dir) )
      maxd += select_db_file( d->d_name, field, chunk );
   if ( !maxd ) {
      closedir( dir );
      return errcode = errMissing;
   }
   db = new SodLcDb [maxd];
   if ( !db ) {
      closedir( dir );
      return errcode = errMemory;
   }

   // open all segments, count all observations, and initialize primary

   primary = 0;
   int o,maxo = 0;
   rewinddir( dir );
   while ( ndb < maxd  &&  (d = readdir(dir)) )
      if ( select_db_file( d->d_name, field, chunk ) ) {
	 sprintf( scratch, "%s/Lightcurve/f_%d/%s", path, field, d->d_name );
	 db[ndb].open( scratch );
	 if ( db[ndb].errcode ) continue;
	 int wjst = db[ndb].hdr->iwjsver < 8 ? 1 : 2;
	 if ( wjstype != wjst )
	    if ( wjstype == -1 ) wjstype = wjst;
	    else continue;
	 if ( db[ndb].hdr->nstars > db[primary].hdr->nstars ) primary = ndb;
	 maxo += db[ndb].hdr->nobs;
	 ++ndb;
      }
   closedir( dir );
   if ( setdb() ) return errcode = db[0].errcode;

   // allocate, fill, and sort the table of contents, then remove duplicates

   if ( !maxo ) return 0;
   toc = new DbTOC [maxo];
   if ( !toc ) return errcode = errMemory;
   for ( int s = 0; s < ndb; ++s )
      for ( o = 0; nobs < maxo  &&  o < db[s].hdr->nobs; ++o )
	 if ( db[s].obs[o].time ) {
	    SodLcDb::Observation &ob = db[s].obs[o];
	    toc[nobs].obsid = ob.obsno;
	    toc[nobs].rtime = ob.isodtime;
	    toc[nobs].segment = s;
	    toc[nobs].obsindex = o;
	    ++nobs;
	 }
   if ( !nobs ) return 0;
   qsort( (char*) toc, nobs, sizeof(SodLcMetaDb::DbTOC)
	  , (int(*)(const void*,const void*)) compare );
   maxo = nobs;
   nobs = 0;
   for ( o = 1; o < maxo; ++o )
      if ( toc[o].obsid != toc[nobs].obsid ) toc[++nobs] = toc[o];
   ++nobs;
   return 0;
}

// The table of contents is to be sorted in ascending order of obsid
// and descending order of reduction time for duplicate observations.
// Duplicate observations are then removed, leaving only the most recent.

int	SodLcMetaDb::compare( const void *i, const void *j )
{
   const DbTOC *I = (const DbTOC*) i;
   const DbTOC *J = (const DbTOC*) j;
   return I->obsid == J->obsid ? J->rtime - I->rtime : I->obsid - J->obsid;
}

int	SodLcMetaDb::obscount( void ) const
{
   return nobs;
}

int	SodLcMetaDb::starcount( void ) const
{
   return ndb ? db[primary].hdr->nstars : 0;
}

SodLcDb::Lightcurve*	SodLcMetaDb::operator[] ( int starindex ) const
{
   if ( segment < 0  ||  segment >= ndb ) return 0;
   return db[segment][starindex];
}

SodLcDb::Lightcurve*	SodLcMetaDb::find( int t, int s ) const
{
   if ( segment < 0  ||  segment >= ndb ) return 0;
   return db[segment].find( t, s );
}

SodLcDb::Lightcurve*	SodLcMetaDb::find( int t, int s, int i, int c ) const
{
   if ( segment < 0  ||  segment >= ndb ) return 0;
   return db[segment].find( t, s, i, c );
}

int	SodLcMetaDb::setdb( int seg )
{
   if ( seg == -1 ) seg = primary;
   if ( seg < 0  ||  seg >= ndb ) return 1;

   segment = seg;
   fname = db[seg].fname;
   hdr = db[seg].hdr;
   obs = db[seg].obs;

   star = 0;
   obsindex = -1;
   index = nobs;

   return 0;
}

int	SodLcMetaDb::next()
{
   if ( index+1 >= nobs ) return obsindex = -1;

   if ( toc[++index].segment != segment )
      for ( star = 0; index < nobs; ++index )
	 if ( toc[index].segment != segment ) {
	    int temp = index;
	    setdb( toc[index].segment );
	    index = temp;
	    if ( star = find( tile, seq ) ) break;
	 }

   return obsindex = star ? toc[index].obsindex : -1;
}

int	SodLcMetaDb::reset( int TILE, int SEQ )
{
   tile = TILE;
   seq = SEQ;
   segment = index = -1;
   star = 0;
   return next();
}

int	SodLcMetaDb::reset( int starindex )
{
   star = 0;
   index = nobs;
   if ( ndb == 0 ) return obsindex = -1;
   SodLcDb::Lightcurve *lc = db[primary][starindex];
   if ( !lc ) return obsindex = -1;
   tile = lc->tile;
   seq = lc->seq;
   segment = index = -1;
   return next();
}
