/* written by Mark Pratt June 1995 */
/* modified by John Doug Reynolds, February 1996 */

#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <malloc.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include "libsite.H"

void LcStat::initialize( void ) {
   map = 0;
   fname = 0;
   clear( );
}   

void LcStat::clear( void ) {
   if( map )
      munmap( map, size );
   map = 0;
   size = 0;
   if( fname )
      free( fname );
   errcode = 0;
   fname = 0;
   Nseq = 0;
}

int LcStat::open( char *path, int tile ) {		// open for reading
   int fd;
   struct stat sbuf;
   char fi[512];

   if( !path )
      return errParams;

   clear( );
   sprintf( fi, "%s/LcSummary/T_%d/Stat_%d", path, tile/100, tile );
   
   if( stat( fi, &sbuf) )
      return errcode = -errno;
   
   if( ( fd = ::open(fi, O_RDONLY) ) == -1 )
      return errcode = -errno;

   size = (int) sbuf.st_size;
   map = (char*) mmap( 0, size, PROT_READ, MAP_PRIVATE, fd, 0);
   close(fd);

   Nseq = size/sizeof( Rec ) - 1;
   fname = strdup( fi );
   _tile = tile;

   if ( size % sizeof( Rec ) )			// funny sized file
      return errcode = errFormat;

   return errNone;
}

int LcStat::open( char *path, int tile, int nrec) { // open to write

   int fd;
   struct stat sbuf;
   char fi[512];

   if( !path )
      return errParams;

   umask( 07 );
   clear( );
   sprintf( fi, "%s/LcSummary/T_%d", path, tile/100 );
   if( stat( fi, &sbuf ) )			// is there a prob with dir?
      if( errno == ENOENT ) {			// is it missing?
	 if( mkdir( fi, 0777 ) )		// make make directory
	    return errcode = -errno;
      } else					// what the hell?
	 return errcode = -errno;

   sprintf( fi, "%s/LcSummary/T_%d/Stat_%d", path, tile/100, tile );
   
   if( ( fd = ::open( fi, O_RDWR | O_CREAT, 0666 ) ) == -1 ) 
      return errcode = -errno;

   if( fstat( fd, &sbuf) )
      return errcode = -errno;

   size = (int) (sbuf.st_size > (nrec+1) * sizeof( Rec ) ? sbuf.st_size :
		 (nrec+1) * sizeof( Rec ) );

   if( ftruncate( fd, size ) )	// stretch to max(size,nrec*recsize)
      return errcode = -errno;
   
   map = (char*) mmap( 0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
   close(fd);
   
   Nseq = size/sizeof( Rec ) - 1;
   fname = strdup( fi );
   _tile = tile;

   return errNone;
}

LcStat::Rec * LcStat::operator[] (int sequence ) const {
   if( !map ) return 0;
   return sequence >= 0 ? sequence <= Nseq ? (Rec *) map + sequence : 0 : 0;
}

int LcStat::print( int sequence, FILE *fp ) {
   if( !(*this)[sequence] )
      return errAccess;

   Rec &rec = *( (*this)[sequence] );
   return print( rec, fp );
}

int LcStat::print( Rec &rec, FILE *fp ) {

   fprintf( fp, "MACHO object %d_%d   Field: %-3d  Chunk: %-3d\n",
	    rec.tile, rec.seq, rec.lot/1000, rec.lot%1000 );
   fprintf( fp, "RA: %s   Dec: %s   cat: %-6d\n",
	    macho_ra_string(rec.RA), macho_dec_string(rec.Dec), rec.cat );
   fprintf( fp, "npoints: %-9d  version: %-8d  statno: %-6d  varno: %-6d\n",
	    rec.npoints, rec.aversion, rec.KGstatnum, rec.KGvarnum );	    
   fprintf( fp, "%39s%30s\n", "RED", "BLUE" );
   fprintf( fp, "%18s: %9d %9.5g %9.5g  %9d %9.5g %9.5g\n",
	    "N mag errave",
	    rec.red.N, rec.red.mag, rec.red.errave,
	    rec.blue.N, rec.blue.mag, rec.blue.errave );
   fprintf( fp, "%18s: %9d %9.5g %9.5g  %9d %9.5g %9.5g\n",
	    "Nr chi2w chi2r",
	    rec.red.Nr, rec.red.chi2w, rec.red.chi2r,
	    rec.blue.Nr, rec.blue.chi2w, rec.blue.chi2r );
   fprintf( fp, "%18s: %9.5g %9.5g %9.5g  %9.5g %9.5g %9.5g\n",
	    "erraver sig magu",
	    rec.red.erraver, rec.red.sig, rec.red.magaveu,
	    rec.blue.erraver, rec.blue.sig, rec.blue.magaveu );
   fprintf( fp, "%18s: %9.5g %9.5g %9.5g  %9.5g %9.5g %9.5g\n",
	    "autocor win[0-1]",
	    rec.red.autoc, rec.red.win[0], rec.red.win[1],
	    rec.blue.autoc, rec.blue.win[0], rec.blue.win[1] );
   fprintf( fp, "%18s: %9.5g %9.5g %9.5g  %9.5g %9.5g %9.5g\n",
	    "\" win[2-4]",
	    rec.red.win[2], rec.red.win[3], rec.red.win[4], 
	    rec.blue.win[2], rec.blue.win[3], rec.blue.win[4] );
   fprintf( fp, "%18s: %9.5g %9.5g %9.5g  %9.5g %9.5g %9.5g\n",
	    "\" win[5-7]",
	    rec.red.win[5], rec.red.win[6], rec.red.win[7], 
	    rec.blue.win[5], rec.blue.win[6], rec.blue.win[7] );
   fprintf( fp, "%18s: %9.5g %9.5g %9.5g  %9.5g %9.5g %9.5g\n",
	    "\" win[8-9] qw",
	    rec.red.win[8], rec.red.win[9], rec.red.qw, 
	    rec.blue.win[8], rec.blue.win[9], rec.blue.qw );
   fprintf( fp, "%18s: %9.5g %9d %9d  %9.5g %9d %9d\n",
	    "qr s5 s5d",
	    rec.red.qr, rec.red.s5, rec.red.s5d,
	    rec.blue.qr, rec.blue.s5, rec.blue.s5d );
   fprintf( fp, "%18s: %9.5g %9.5g %9.5g  %9.5g %9.5g %9.5g\n",
	    "decor_c, s miss",
	    rec.red.decorcon, rec.red.decorslope, rec.red.missave,
	    rec.blue.decorcon, rec.blue.decorslope, rec.blue.missave );
   fprintf( fp, "%18s: %9.5g %9.5g %9.5g  %9.5g %9.5g %9.5g\n",
	    "dsky dskysig crd",
	    rec.red.dskyave, rec.red.dskysig, rec.red.crdave,
	    rec.blue.dskyave, rec.blue.dskysig, rec.blue.crdave );
   fprintf( fp, "%18s: %9.5g %9.5g %9.5g  %9.5g %9.5g %9.5g\n",
	    "psf psfsig cosm",
	    rec.red.psfave, rec.red.psfsig, rec.red.cosmave,
	    rec.blue.psfave, rec.blue.psfsig, rec.blue.cosmave );
   fprintf( fp, "%18s: %9d %9.5g %9.5g  %9d %9.5g %9.5g\n",
	    "ml_0 nhip max time",
	    rec.red.m1stat[0].nhipoints, rec.red.m1stat[0].filmax,
	    rec.red.m1stat[0].filtime,  rec.blue.m1stat[0].nhipoints,
	    rec.blue.m1stat[0].filmax, rec.blue.m1stat[0].filtime );
   fprintf( fp, "%18s: %9d %9.5g %9.5g  %9d %9.5g %9.5g\n",
	    "ml_1 nhip max time",
	    rec.red.m1stat[1].nhipoints, rec.red.m1stat[1].filmax,
	    rec.red.m1stat[1].filtime,  rec.blue.m1stat[1].nhipoints,
	    rec.blue.m1stat[1].filmax, rec.blue.m1stat[1].filtime );
   fprintf( fp, "%18s: %9d %9.5g %9.5g  %9d %9.5g %9.5g\n",
	    "ml_2 nhip max time",
	    rec.red.m1stat[2].nhipoints, rec.red.m1stat[2].filmax,
	    rec.red.m1stat[2].filtime,  rec.blue.m1stat[2].nhipoints,
	    rec.blue.m1stat[2].filmax, rec.blue.m1stat[2].filtime );
   fprintf( fp, "\n" );
   fprintf( fp, "%18s: %9d %9.5g %9.5g %9.5g\n",
	    "N B-R B-Rsig chi2w",
	    rec.br.N, rec.br.ave, rec.br.sig, rec.br.chi2w );
   fprintf( fp, "%18s: %11.4g\n", "probspear" , rec.br.probspear );
   fprintf( fp, "%18s: %11.4g\n", "cross-corr", rec.br.cross );

   return SiteErr::errNone;
}
