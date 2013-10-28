/* written by John Doug Reynolds, February 1998 */

#include "libsite.H"
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#ifndef _REENTRANT
#define strtok_r(s,d,n) strtok(s,d)
#endif

/*	MACHO BadObs file format specification

	A comment line contains either only whitespace characters or #
	as the first non-whitespace character.  All other lines must
	be well formed badobs descriptions.  A badobs description is a
	range followed by an optional amp list, and optionally
	concluded with a comment introduced by #.

	A range is either obsid or [obsid1]-[obsid2].  In the second
	form obsid1 defaults to the smallest integer and obsid2
	defaults to the largest integer, and the range thus specified
	includes both endpoints.  obsid2 may be less than obsid1.
	Negative observation numbers are not supported.

	The amp list defaults to all amps.  When specified, the list
	is a single token consisting of any set of hex digits in upper
	or lower case.  Order is not important and duplicates are
	tolerated.  Each digit refers to a single bad amplifier.

	The order of badobs descriptions is not important.  In other
	words, the BadObs file is expected to contain ranges specified
	in random order and overlapping each other in arbitrary ways.

	#Example:
	300			# single bad observation
	50000-			# remove all observations from 50000 on
	40000-41000	EF	# temporary ccd 7 failure
*/

MachoBadObs::MachoBadObs()
{
   count = maxcap = 0;
   list = 0;
}

MachoBadObs::~MachoBadObs()
{
   if ( list ) delete [] list;
}

int	MachoBadObs::capacity( void ) const
{
   return maxcap;
}

void	MachoBadObs::extend( int cap )
{
   if ( cap > maxcap ) {
      struct Range *newlist = new Range [cap];
      if ( list ) {
	 memcpy( newlist, list, count * sizeof(Range) );
	 delete [] list;
      }
      list = newlist;
      maxcap = cap;
   }
}

static int	rangecompare( const void *a_, const void *b_ )
{
   int a = ((const MachoBadObs::Range*)a_)->min;
   int b = ((const MachoBadObs::Range*)b_)->min;
   return a >= b ? a == b ? 0 : 1 : -1;
}

void	MachoBadObs::normalize( void )
{
   if ( count <= 0 ) return;

   // sort list by Range::min and save as local variable, reset list

   qsort( list, count, sizeof(Range), rangecompare );

   int i, j, nrange = count;
   struct Range *range = list;

   list = 0;
   count = maxcap = 0;
   extend( nrange );

   // rebuild list as sequence of non-overlapping ranges

   for ( i = 0; i < nrange; ) {
      struct Range r = range[i];

      // search for minimum upper boundary of current range

      for ( j = i+1; j < nrange; ++j )
	 if ( range[j].min == r.min ) {
	    if ( range[j].max < r.max ) r.max = range[j].max;
	    r.amps |= range[j].amps;
	 }
	 else {
	    if ( range[j].min <= r.max ) r.max = range[j].min - 1;
	    break;
	 }

      // truncate contributing ranges and delete any that become empty

      if ( r.max < INT_MAX )
	 for ( j = i; j < nrange && range[j].min == r.min; ++j ) {
	    if ( (range[j].min = r.max + 1) > range[j].max )
	       range[j] = range[i++];
	 }
      else i = nrange;

      // extend the previous range or start a new one

      if ( r.amps )
	 if ( count > 0 &&
	      list[count-1].max + 1 == r.min && list[count-1].amps == r.amps )
	    list[count-1].max = r.max;
	 else {
	    if ( count >= maxcap ) extend( count < 256 ? 256 : 2 * count );
	    list[count++] = r;
	 }
   }

   delete [] range;
}

int	MachoBadObs::load( const char *filename, int verbose )
{
   count = 0;

   int line = 0;
   int error = 0;
   char buf[1024];
   FILE *fi = fopen( filename, "r" );
   if ( ! fi ) {
      if ( verbose ) perror(filename);
      return -errno;
   }

   while ( !error && fgets(buf,1024,fi) ) {
      ++line;

      const char WHITESPACE[] = " \t\n";
      char *next, *tok = strtok_r(buf,WHITESPACE,&next);
      if ( !tok || *tok == '#' ) continue;

      // parse the observation range token

      struct Range r;
      int n, len = strlen(tok);

      if ( sscanf(tok,"%d-%d%n",&r.min,&r.max,&n) != 2 ) {
	 if ( *tok == '-' ) {
	    r.min = INT_MIN;
	    if ( len == 1 ) r.max = INT_MAX;
	    else error = sscanf(tok,"-%d%n",&r.max,&n) != 1 || n != len;
	 }
	 else if ( sscanf(tok,"%d%n",&r.min,&n) == 1 )
	    if ( n == len ) r.max = r.min;
	    else if ( tok[n] == '-' && n+1 == len ) r.max = INT_MAX;
	    else error = 1;
	 else
	    error = 1;
      }
      else {
	 if ( r.max < r.min ) { int t = r.min; r.min = r.max; r.max = t; }
	 error = n != len;
      }

      if ( error ) continue;

      // parse the amplifier list

      tok = strtok_r(NULL,WHITESPACE,&next);
      if ( tok && *tok != '#' ) {
	 r.amps = 0;
	 while ( !error && *tok )
	    if ( sscanf(tok++,"%1x",&n) == 1 ) r.amps += 1 << n;
	    else error = 1;
	 if ( ! error ) {
	    tok = strtok_r(NULL,WHITESPACE,&next);
	    error = tok && *tok != '#';
	 }
      }
      else r.amps = 0xFFFF;

      // append the new range to list

      if ( ! error ) {
	 if ( count >= maxcap ) extend( count < 256 ? 256 : 2 * count );
	 list[count++] = r;
      }
   }

   fclose( fi );

   if ( error ) {
      count = 0;
      if ( verbose )
	 fprintf( stderr, "%s: format error at line %d\n", filename, line );
      return SiteErr::errFormat;
   }

   normalize();
   return 0;
}

int	MachoBadObs::find( int obsid ) const
{
   int i0 = 0;
   int i1 = count - 1;

   while ( i1 >= i0 ) {
      int i = (i0 + i1) / 2;
      if ( obsid < list[i].min ) i1 = i - 1;
      else if ( obsid > list[i].max ) i0 = i + 1;
      else return i;
   }

   return -1;
}

int	MachoBadObs::amps( int obsid, int &index ) const
{
   if ( count <= 0 ) return 0;

   if ( index < 0 )
      index = 0;
   else if ( index >= count )
      index = count - 1;

   if ( obsid < list[index].min ) {
      while ( index > 0 && obsid <= list[index-1].max )
	 if ( obsid >= list[--index].min ) return list[index].amps;
      return 0;
   }

   if ( obsid > list[index].max ) {
      while ( index+1 < count && obsid >= list[index+1].min )
	 if ( obsid <= list[++index].max ) return list[index].amps;
      return 0;
   }

   return list[index].amps;
}

int	MachoBadObs::amps( int obsid ) const
{
   int index = find(obsid);
   return index == -1 ? 0 : list[index].amps;
}
