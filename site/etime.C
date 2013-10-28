/* written by John Doug Reynolds, March 1997 */

/*	ExplicitTM	an extension of struct tm
 *
 * The purpose of this class is to provide a sane set of tools for
 * handling the ASCII representation of a time when the time zone is
 * explicit.  The underlying system calls handle this situation with
 * unusually poor grace.  Several of these methods allow an optional
 * pointer to a time zone name.  In each case the default is a null
 * pointer, and the value of the environment variable TZ is used.  If
 * TZ is not in the environment, or a pointer to string of length zero
 * is provided, the local default time zone is used.  If TZ is in the
 * environment, but has no value, GMT is used.
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "libsite.H"

#if defined __sunos4__ || defined __osf1__
#define __bsdtime__
#endif

#ifdef __linux__
/*extern long int altzone;*/
#endif

#ifdef __sunos4__
extern "C" time_t mktime( struct tm * );
#endif

static char *temptz = 0;	// scratch space for the environment

static const char *meridian[] = {
   "AM", "PM"
};
static const char *weekday[] = {
   "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"
};
static const char *weekday3[] = {
   "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"
};
static const char *month[] = {
   "January", "February", "March", "April", "May", "June",
   "July", "August", "September", "October", "November", "December"
};
static const char *month3[] = {
   "Jan", "Feb", "Mar", "Apr", "May", "Jun",
   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};

ExplicitTM::ExplicitTM( void )
{
   out = 0;
   maxout = 0;
   if ( !temptz ) {
      temptz = new char [ZoneMax + 3];
      putenv( strcpy( temptz, "_tz=" ) );
   }
   reset("");
}

ExplicitTM::~ExplicitTM( void )
{
   if ( out ) delete [] out;
}

const char*	ExplicitTM::etime( const time_t &clock, const char *tzone )
{
   return reset(clock,tzone) ? 0 : format( "%a %b %d %X %Z %Y" );
}

int	ExplicitTM::reset( const char *tzone )
{
   memset( this, 0, sizeof(struct tm) );
   if ( out ) out[0] = '\0';
   zone3[0] = zone[0] = '\0';
   gmtoff = 0;

   if ( !tzone ) {
      tzone = getenv( "TZ" );
      if ( tzone && !*tzone ) tzone = "GMT";
   }
   if ( tzone ) {
      if ( strlen( tzone ) >= ZoneMax ) return 1;
      strcpy( zone, tzone );
   }

   return 0;
}

int	ExplicitTM::reset( const time_t &clock, const char *tzone )
{
   extern char **environ;

   if ( reset( tzone ) ) return 1;

   // find and obscure the original timezone
   int tz;
   for ( tz = 0; environ[tz]; ++tz )
      if ( strncmp(environ[tz],"TZ=",3) == 0 ) break;
   if ( environ[tz] ) environ[tz][0] = '\0';

   // convert the time using the temporary timezone
   if ( zone[0] ) sprintf( temptz, "TZ=%s", zone );
   tzset();
   struct tm *tmp = localtime( &clock );
   memcpy( this, tmp, sizeof(struct tm) );
#ifdef __bsdtime__
   char *z3 = tm_zone;
   gmtoff = (int) tm_gmtoff;
#else
   char *z3 = (tm_isdst == 0 || tm_isdst == 1) ? tzname[ tm_isdst ] : 0;
   /*gmtoff = (tm_isdst == 1) ? -altzone : -timezone;*/
   gmtoff = -timezone;
#endif
   if ( !gmtoff ) {
      strcpy( zone, "GMT" );
      memcpy( zone3, zone, 3 );
   }
   else if ( z3 ) memcpy( zone3, z3, 3 );

   // restore the original timezone
   if ( zone[0] ) strcpy( temptz, "_tz=" );
   if ( environ[tz] ) environ[tz][0] = 'T';
   tzset();

   return 0;
}

int	ExplicitTM::get( time_t &clock )
{
   extern char **environ;

   // find and obscure the original timezone
   int tz;
   for ( tz = 0; environ[tz]; ++tz )
      if ( strncmp(environ[tz],"TZ=",3) == 0 ) break;
   if ( environ[tz] ) environ[tz][0] = '\0';

   // convert the time using the temporary timezone
   if ( zone[0] ) sprintf( temptz, "TZ=%s", zone );
   tzset();
   tm_isdst = -1;
   clock = mktime( this );
#ifdef __bsdtime__
   char *z3 = tm_zone;
   gmtoff = (int) tm_gmtoff;
#else
   char *z3 = (tm_isdst == 0 || tm_isdst == 1) ? tzname[ tm_isdst ] : 0;
   /*gmtoff = (tm_isdst == 1) ? -altzone : -timezone;*/
   gmtoff = -timezone;
#endif
   memcpy( zone3, z3 ? z3 : "\0\0\0", 3 );

   // restore the original timezone
   if ( zone[0] ) strcpy( temptz, "_tz=" );
   if ( environ[tz] ) environ[tz][0] = 'T';
   tzset();

   return clock == -1 || (!gmtoff && strncmp(zone3,"GMT",3));
}

int	ExplicitTM::format( const char *fmt, int &len )
{
   char scratch[32];

   while ( *fmt ) {
      int index, nbytes = -1;
      const char *buf = scratch;

      if ( *fmt == '%' ) {
	 int code = *++fmt;
	 if ( code ) ++fmt;
	 else continue;

	 switch ( code ) {
	 case 'a':
	 case 'A':
	    index = tm_wday % 7;
	    if ( index < 0 ) index += 7;
	    buf = code == 'a' ? weekday3[index] : weekday[index];
	    break;
	 case 'b':
	 case 'B':
	    index = tm_mon % 12;
	    if ( index < 0 ) index += 12;
	    buf = code == 'b' ? month3[index] : month[index];
	    break;
	 case 'c':
	    if ( format("%x %X",len) ) return 1;
	    continue;
	 case 'd':
	    sprintf( scratch, "%02d", tm_mday );
	    break;
	 case 'G':
	    sprintf( scratch, "GMT%+.2g", gmtoff / 3600. );
	    break;
	 case 'H':
	    sprintf( scratch, "%02d", tm_hour );
	    break;
	 case 'I':
	    index = tm_hour % 12;
	    if ( !index ) index = 12;
	    sprintf( scratch, "%02d", index );
	    break;
	 case 'j':
	    sprintf( scratch, "%03d", tm_yday + 1 );
	    break;
	 case 'm':
	    sprintf( scratch, "%2d", tm_mon + 1 );
	    break;
	 case 'M':
	    sprintf( scratch, "%02d", tm_min );
	    break;
	 case 'p':
	    index = tm_hour % 24;
	    if ( index < 0 ) index += 24;
	    buf = meridian[ index / 12 ];
	    break;
	 case 'S':
	    sprintf( scratch, "%02d", tm_sec );
	    break;
	 case 'w':
	    sprintf( scratch, "%d", tm_wday );
	    break;
	 case 'x':
	    if ( format("%m/%d/%y",len) ) return 1;
	    break;
	 case 'X':
	    if ( format("%H:%M:%S",len) ) return 1;
	    continue;
	 case 'y':
	    sprintf( scratch, "%02d", tm_year % 100 );
	    break;
	 case 'Y':
	    sprintf( scratch, "%d", tm_year + 1900 );
	    break;
	 case 'z':
	    if ( zone3[0] ) {
	       buf = zone3;
	       nbytes = 3;
	    }
	    else continue;
	    break;
	 case 'Z':
	    if ( zone[0] ) buf = zone;
	    else if ( format("%G",len) ) return 1;
	    else continue;
	    break;
	 case '%':
	    scratch[0] = code;
	    nbytes = 1;
	    break;
	 default:
	    scratch[0] = '%';
	    scratch[1] = code;
	    nbytes = 2;
	    break;
	 }
      }
      else {
	 const char *end = strchr( fmt, '%' );
	 if ( end )
	    nbytes = end - fmt;
	 else
	    nbytes = strlen( fmt );
	 buf = fmt;
	 fmt += nbytes;
      }

      if ( nbytes < 0 ) nbytes = strlen( buf );
      if ( len + nbytes >= maxout ) {
	 char *temp = new char [maxout += 256];
	 if ( out ) {
	    memcpy( temp, out, len );
	    delete [] out;
	 }
	 out = temp;
      }
      if ( nbytes ) {
	 memcpy( out+len, buf, nbytes );
	 len += nbytes;
      }
      if ( out ) out[len] = '\0';
   }

   return 0;
}

const char*	ExplicitTM::format( const char *fmt )
{
   if ( !fmt ) return out;

   int len = 0;
   if ( out ) out[0] = '\0';
   return format(fmt,len) ? 0 : out;
}

static int parse( const char* &str, const char *dict[], int max )
{
   int i, n;
   while ( isspace(*str) ) ++str;
   for ( i = 0; i < max; ++i ) {
      n = strlen( dict[i] );
      if ( strncmp( str, dict[i], n ) == 0 ) break;
   }
   if ( i < max ) {
      str += n;
      return i;
   }
   else return -1;
}

int	ExplicitTM::scan( const char* &str, const char *fmt )
{
   char scratch[16];
   int index, nbytes;

   while ( *fmt )
      if ( *fmt == '%' ) {
	 int code = *++fmt;
	 if ( code ) ++fmt;
	 else continue;

	 switch ( code ) {
	 case 'a':
	 case 'A':
	    index = parse( str, code == 'a' ? weekday3 : weekday, 7 );
	    if ( index < 0 ) return 1;
	    tm_wday = index;
	    break;
	 case 'b':
	 case 'B':
	    index = parse( str, code == 'b' ? month3 : month, 12 );
	    if ( index < 0 ) return 1;
	    tm_mon = index;
	    break;
	 case 'c':
	    if ( scan( str, "%x %X" ) ) return 1;
	    break;
	 case 'd':
	    if ( sscanf( str, "%d%n", &tm_mday, &nbytes ) != 1 ) return 1;
	    str += nbytes;
	    break;
	 case 'H':
	 case 'I':
	    if ( sscanf( str, "%d%n", &tm_hour, &nbytes ) != 1 ) return 1;
	    if ( code == 'I' && tm_hour == 12 ) tm_hour = 0;
	    str += nbytes;
	    break;
	 case 'j':
	    if ( sscanf( str, "%d%n", &tm_yday, &nbytes ) != 1 ) return 1;
	    str += nbytes;
	    --tm_yday;
	    break;
	 case 'm':
	    if ( sscanf( str, "%d%n", &tm_mon, &nbytes ) != 1 ) return 1;
	    str += nbytes;
	    --tm_mon;
	    break;
	 case 'M':
	    if ( sscanf( str, "%d%n", &tm_min, &nbytes ) != 1 ) return 1;
	    str += nbytes;
	    break;
	 case 'p':
	    if ( (index = parse( str, meridian, 2 )) < 0 ) return 1;
	    if ( index == 1 ) tm_hour += 12;
	    break;
	 case 'S':
	    if ( sscanf( str, "%d%n", &tm_sec, &nbytes ) != 1 ) return 1;
	    str += nbytes;
	    break;
	 case 'w':
	    if ( sscanf( str, "%d%n", &tm_wday, &nbytes ) != 1 ) return 1;
	    str += nbytes;
	    break;
	 case 'x':
	    if ( scan( str, "%m/%d/%y" ) ) return 1;
	    break;
	 case 'X':
	    if ( scan( str, "%H:%M:%S" ) ) return 1;
	    break;
	 case 'y':
	    if ( sscanf( str, "%d%n", &tm_year, &nbytes ) != 1 ) return 1;
	    if ( tm_year < 20 ) tm_year += 100;
	    str += nbytes;
	    break;
	 case 'Y':
	    if ( sscanf( str, "%d%n", &tm_year, &nbytes ) != 1 ) return 1;
	    tm_year -= 1900;
	    str += nbytes;
	    break;
	 case 'z':
	 case 'Z':
	 case 'G':
	    sprintf( scratch, "%%%ds%%n", ZoneMax-1 );
	    if ( sscanf( str, scratch, zone, &nbytes ) != 1 ) return 1;
	    str += nbytes;
	    break;
	 case '%':
	    if ( *str != code ) return 1;
	    ++str;
	    break;
	 default:
	    return 1;
	 }
      }
      else {
	 const char *end = strchr( fmt, '%' );
	 if ( end )
	    nbytes = end - fmt;
	 else
	    nbytes = strlen( fmt );
	 if ( strncmp( str, fmt, nbytes ) ) return 1;
	 str += nbytes;
	 fmt += nbytes;
      }

   return 0;
}

const char*	ExplicitTM::sscan( const char* str, const char *fmt )
{
   return ( !str || !fmt || scan(str,fmt) ) ? 0 : str;
}
