/* written by John Doug Reynolds, March 1998 */

#include <ctype.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include "libsite.H"

int usage( const char *appname )
{
   fprintf( stderr, "\
usage %s: [-v] [-q <count>] [-m <count>] [-u <count>] [-f <count>]\n\
              [-min <obsid>] [-max <obsid>] [-c <field>.<tile> <list>]\n\
              [-x <sodset> <list>] [-r <sodset>] [-t <contents>]\n\
              [-p [<coords>] <lc>] <filename>\n\n\
   Creates and manipulates portable lightcurve databases.\n\n\
   These options are meaningful only in combination with -c below.  They\n\
   control access to archived photometry databases.\n\n\
     -q <count>     Enable archival access for only the latest <count>\n\
                    observations.  Default is zero, meaning never disable.\n\
     -m <count>     Allow no more than <count> database fragment merges.\n\
                    Unlimited by default.\n\
     -u <count>     Allow no more than <count> database uncompressions.\n\
                    Unlimited by default.\n\
     -f <count>     Allow no more than <count> nearline database retrievals.\n\
                    Unlimited by default.\n\
     -s             Skip inaccessible databases and allow quota violations.\n\
     -v             Print diagnostic messages from the extractor.\n\n\
   These options apply more generally.\n\n\
     -min <obsid>   Exclude observations before <obsid>.  Default is zero.\n\
     -max <obsid>   Exclude observations after <obsid>.  Default is zero,\n\
                    meaning the latest observation.\n\n\
   These options select the mode of operation, and are mutually exclusive.\n\
   The default mode just prints the header of the specified database.\n\n\
     -c <field>.<tile> <list>\n\
                    Create a new database containing the indicated stars.\n\
     -x <sodset> <list>\n\
                    Create an excerpt of <sodset> with the indicated stars.\n\
     -r <sodset>\n\
                    Create a copy of <sodset> in reverse byte-order.\n\
     -t stars | obsids | chunks\n\
                    Print the indicated table of contents.\n\
     -p [chip | amp | fits | none] <sequence> | %%<index>\n\
                    Print the lightcurve indicated by sequence number or\n\
                    table-of-contents index, such as %%0.  Pixel coordinates\n\
                    are represented in the amp system by default.\n\n\
   <list> refers to a list of stars in one of the following formats:\n\n\
        -           Read an unordered list of starids from standard input.\n\
                    One <field>.<tile>.<seqn> per line.  Comments ignored.\n\
                    Duplicates and field or tile mismatches are removed.\n\
        <seqn>[,<seqn>...]\n\
                    Unordered list of sequence numbers.  Duplicates removed.\n\
        <skip>/<count>[/<inc>]\n\
                    Take every <inc>th star, to a maximum of <count> stars,\n\
                    starting after the <skip>th star.  <inc> defaults to 1.\n\
        all         Equivalent to 0/65536.\n\n\
"
	    , appname );
   return 1;
}

int open( class SodophotSet &ss, const char *filename )
{
   if ( ss.open(filename) ) {
      if ( ss.version() < 0 )
	 fprintf( stderr, "%s might be readable once byte-swapped\n"
		  , filename );
      else ss.perror( filename );
      return ss.errcode;
   }
   else return 0;
}

int header( const char *filename )
{
   SodophotSet ss;
   if ( open(ss,filename) ) return ss.errcode;

   int i;
   SodophotSet::Header &h = *ss.hdr;
   printf( "sodset_%d.%d_%c\n", h.field, h.tile, h.source );
   printf( "%20s%12s%12s%12s%12s\n"
	   , "", "contents", "capacity", "min id", "max id" );
   printf( "   format %d\n", ss.version() );
   printf( "%20s%12d%12d", "stars", h.count.s, ss.maximum().s );
   if ( h.count.s )
      printf( "%12d%12d", ss.star[0].seqn, ss.star[h.count.s-1].seqn );
   printf( "\n" );
   printf( "%20s%12d%12d", "observations", h.count.o, ss.maximum().o );
   if ( h.count.o )
      printf( "%12d%12d", ss.obs[0].obsid, ss.obs[h.count.o-1].obsid );
   printf( "\n" );
   printf( "%20s%12d%12d", "chunks", h.count.c, ss.maximum().c );
   for ( i = 0; i < 128 && !ss.chnk[i]; ++i );
   if ( i < 128 ) printf( "%12d", i );
   for ( i = 127; i >= 0 && !ss.chnk[i]; --i );
   if ( i >= 0 ) printf( "%12d", i );
   printf( "\n\n" );

   return 0;
}

int ushortcompare( const void *i, const void *j )
{
   return (int) (*((const unsigned short*)i) - *((const unsigned short*)j));
}

int getlist( const char *str, int field, int tile
	     , unsigned short *&list, int &skip, int &count, int &inc )
{
   list = 0; 
   skip = count = 0;
   inc = 1;

   // check for implicit list specification

   if ( strcmp( str, "all" ) == 0 ) {
      count = 65536;
      return 0;
   }
   else if ( strchr( str, '/' ) ) {
      if ( sscanf( str, "%d/%d/%d", &skip, &count, &inc ) < 2
	   || skip < 0 || count < 0 || inc < 1 )
	 return SiteErr(SiteErr::errParams).perror(str);
      return 0;
   }

   // construct a list from str or stdin

   int max = 0;
   if ( strcmp( str, "-" ) == 0 ) str = 0;

   for (;;) {
      unsigned short seqn;

      if ( str ) {
	 int n;
	 if ( ! *str ) break;
	 if ( sscanf( str, "%hu%n", &seqn, &n ) != 1 ) {
	    fprintf( stderr, "improper sequence list: %s\n", str );
	    return 1;
	 }
	 if ( *(str+=n) == ',' ) ++str;
      }
      else {
	 int f, t;
	 char buf[1024], id[18];
	 if ( ! fgets( buf, 1024, stdin ) ) break;
	 if ( sscanf( buf, "%17s", id ) != 1 || strchr( "#!", *id ) ) continue;
	 if ( sscanf( id, "%d.%d.%hu", &f, &t, &seqn ) != 3 ) {
	    fprintf( stderr, "improper starid: %s\n", id );
	    return 1;
	 }
	 if ( f != field || t != tile ) continue;
      }

      if ( count >= max ) {
	 unsigned short *temp = new unsigned short [max+512];
	 if ( ! temp ) return SiteErr(SiteErr::errMemory).perror("getlist");
	 if ( list ) {
	    memcpy( temp, list, count * sizeof(unsigned short) );
	    delete [] list;
	 }
	 list = temp;
	 max += 512;
      }

      list[count++] = seqn;
   }

   // sort and remove duplicates

   if ( count ) {
      int i = 0, j = 0, last = -1;
      qsort( list, count, sizeof(unsigned short), ushortcompare );
      for (; i < count; ++i )
	 if ( list[i] != last ) last = list[j++] = list[i];
      count = j;
   }

   return 0;
}

int extract( const char *filename, struct ExtractLcOpts &opts
	     , const char *fieldtile, const char *starlist )
{
   if ( sscanf( fieldtile, "%d.%d", &opts.field, &opts.tile ) != 2 )
      return SiteErr(SiteErr::errParams).perror(fieldtile);

   int skip, count, inc;
   unsigned short *list;
   if ( getlist( starlist, opts.field, opts.tile, list, skip, count, inc ) )
      return 1;

   if ( count == 0 ) {
      fprintf( stderr, "no stars from tile %d.%d specified\n"
	       , opts.field, opts.tile );
      return 1;
   }

   SodophotSet ss;
   if ( list ? ss.extract(opts,count,list) : ss.extract(opts,skip,count,inc) )
      return ss.perror();

   if ( ss.errcode = ss.write( filename ) ) return ss.perror( filename );
   printf( "%s: %d stars by %d observations\n"
	   , filename, ss.hdr->count.s, ss.hdr->count.o );

   return 0;
}

int excerpt( const char *filename, const ExtractLcOpts &opts
	     , const char *sodset, const char *starlist )
{
   SodophotSet ss;
   if ( open(ss,sodset) ) return ss.errcode;

   int *index;
   int skip, count, inc;
   unsigned short *list;
   if ( getlist( starlist, ss.hdr->field, ss.hdr->tile
		 , list, skip, count, inc ) ) return 1;

   if ( list ) {
      int i, j, max = count;
      index = new int [max];
      for ( i = j = count = 0; i < ss.hdr->count.s && j < max; )
	 if ( ss.star[i].seqn >= list[j] )
	    if ( ss.star[i].seqn == list[j++] )
	       index[count++] = i++;
	    else
	       list[j-count-1] = list[j-1];
	 else ++i;
      if ( count < max ) {
	 fprintf( stderr, "following %d stars missing from %s\n"
		  , max - count, sodset );
	 for ( i = 0; i < max - count; ++i )
	    fprintf( stderr, "\t%d.%d.%d\n"
		     , ss.hdr->field, ss.hdr->tile, list[i] );
      }
   }
   else {
      int i = skip, total = ((ss.hdr->count.s - skip) + (inc - 1)) / inc;
      if ( total > count ) total = count;
      index = new int [total];
      for ( count = 0; count < total; i += inc ) index[count++] = i;
   }

   if ( count == 0 ) {
      fprintf( stderr, "no stars specified from %s\n", sodset );
      return 1;
   }

   if ( ss.errcode = ss.write( filename, count, index ) )
      return ss.perror( filename );
   printf( "%s: %d stars by %d observations\n"
	   , filename, count, ss.hdr->count.o );

   return 0;
}

int contents( const char *filename, const ExtractLcOpts &opts, const char *id )
{
   int i;
   SodophotSet ss;
   if ( open(ss,filename) ) return ss.errcode;

   if ( strcmp( id, "stars" ) == 0 )
      for ( i = 0; i < ss.hdr->count.s; i += 13 ) {
	 for ( int j = i; j < i+13 && j < ss.hdr->count.s; ++j )
	    printf( "%6d", ss.star[j].seqn );
	 printf( "\n" );
      }
   else if ( strcmp( id, "obsids" ) == 0 ) {
      int max = opts.maxobs ? opts.maxobs : INT_MAX;
      for ( i = 0; i < ss.hdr->count.o && ss.obs[i].obsid < opts.minobs; ++i );
      for (; i < ss.hdr->count.o && ss.obs[i].obsid <= max; i += 11 ) {
	 for ( int j = i; j < i+11 && j < ss.hdr->count.o; ++j )
	    if ( ss.obs[j].obsid <= max ) printf( "%7d", ss.obs[j].obsid );
	 printf( "\n" );
      }
   }
   else if ( strcmp( id, "chunks" ) == 0 )
      for ( i = 0; i < 128; ++i ) {
	 if ( ss.chnk[i] ) printf( "%d\n", i );
      }

   return 0;
}

int display( const char *filename, const ExtractLcOpts &opts
	     , const char *coords, const char *lc )
{
   int code;
   if ( strcmp( coords, "chip" ) == 0 )
      code = 'c';
   else if ( strcmp( coords, "amp" ) == 0 )
      code = 'a';
   else if ( strcmp( coords, "fits" ) == 0 )
      code = 'f';
   else if ( strcmp( coords, "none" ) == 0 )
      code = 0;
   else
      return SiteErr(SiteErr::errParams).perror(coords);

   SodophotSet ss;
   if ( open(ss,filename) ) return ss.errcode;

   int o = 0;
   while ( o < ss.hdr->count.o && ss.obs[o].obsid < opts.minobs ) ++o;

   int starindex = (*lc == '%') ? atoi(lc+1) : ss.find(atoi(lc));
   if ( starindex < 0 || starindex >= ss.hdr->count.s )
      return SiteErr(SiteErr::errParams).perror(lc);

   int obsindex = o, count = 0;
   if ( opts.maxobs )
      while ( o < ss.hdr->count.o && ss.obs[o++].obsid <= opts.maxobs )
	 ++count;
   else
      count = -1;

   return ss.printlc( stdout, starindex, code, obsindex, count );
}

int main( int argc, char *argv[] )
{
   int i, mode = 0;
   char *arg1, *arg2, *filename = 0;
   struct ExtractLcOpts opts;

   opts.merge = opts.unzip = opts.fetch = INT_MAX;
   opts.active = getenv( "MACHOACTIVE" );
   opts.cache = getenv( "MACHOLCCACHE" );

   for ( i = 1; i < argc; ++i )
      if ( strcmp( argv[i], "-h" ) == 0 )
	 return !usage( argv[0] );
      else if ( strcmp( argv[i], "-v" ) == 0 )
	 opts.verbose = 1;
      else if ( strcmp( argv[i], "-s" ) == 0 )
	 opts.tolerant = 1;
      else if ( strcmp( argv[i], "-q" ) == 0 && i+1 < argc )
	 opts.quota = atoi( argv[++i] );
      else if ( strcmp( argv[i], "-m" ) == 0 && i+1 < argc )
	 opts.merge = atoi( argv[++i] );
      else if ( strcmp( argv[i], "-u" ) == 0 && i+1 < argc )
	 opts.unzip = atoi( argv[++i] );
      else if ( strcmp( argv[i], "-f" ) == 0 && i+1 < argc )
	 opts.fetch = atoi( argv[++i] );
      else if ( strcmp( argv[i], "-min" ) == 0 && i+1 < argc )
	 opts.minobs = atoi( argv[++i] );
      else if ( strcmp( argv[i], "-max" ) == 0 && i+1 < argc )
	 opts.maxobs = atoi( argv[++i] );
      else if ( ! mode && strcmp( argv[i], "-c" ) == 0 && i+2 < argc ) {
	 mode = 'c';
	 arg1 = argv[++i];
	 arg2 = argv[++i];
      }
      else if ( ! mode && strcmp( argv[i], "-x" ) == 0 && i+2 < argc ) {
	 mode = 'x';
	 arg1 = argv[++i];
	 arg2 = argv[++i];
      }
      else if ( ! mode && strcmp( argv[i], "-r" ) == 0 && i+1 < argc ) {
	 mode = 'r';
	 arg1 = argv[++i];
      }
      else if ( ! mode && strcmp( argv[i], "-t" ) == 0 && i+1 < argc ) {
	 mode = 't';
	 arg1 = argv[++i];
      }
      else if ( ! mode && strcmp( argv[i], "-p" ) == 0 && i+1 < argc
		&& (isdigit(argv[i+1][0]) || argv[i+1][0] == '%') ) {
	 mode = 'p';
	 arg1 = "amp";
	 arg2 = argv[++i];
      }
      else if ( ! mode && strcmp( argv[i], "-p" ) == 0 && i+2 < argc ) {
	 mode = 'p';
	 arg1 = argv[++i];
	 arg2 = argv[++i];
      }
      else if ( ! filename )
	 filename = argv[i];
      else
	 return usage( argv[0] );

   if ( ! filename ) return usage( argv[0] );

   switch ( mode ) {
   default: return header( filename );
   case 'c': return extract( filename, opts, arg1, arg2 );
   case 'x': return excerpt( filename, opts, arg1, arg2 );
   case 'r': return sodset_byte_reverse( filename, arg1 );
   case 't': return contents( filename, opts, arg1 );
   case 'p': return display( filename, opts, arg1, arg2 );
   }
}
