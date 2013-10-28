/* written by John Doug Reynolds, October 1997 */

#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "Macho.h"
#include "libsite.H"

int	SiteErr::perror( const char *label ) const
{
   char *estr;

   switch ( errcode ) {
   default:
      if ( !(estr = strerror(-errcode)) ) estr = "unknown error"; break;
   case errNone:
      estr = "no error"; break;
   case errProgram:
      estr = "internal program error"; break;
   case errParams:
      estr = "bad parameters"; break;
   case errEnviron:
      estr = "missing environment variable"; break;
   case errMemory:
      estr = "memory allocation error"; break;
   case errAccess:
      estr = "resource access error"; break;
   case errFormat:
      estr = "resource format error"; break;
   case errMissing:
      estr = "resource not found"; break;
   case errVersion:
      estr = "unimplemented logic"; break;
   case errMath:
      estr = "math exception"; break;
   case errCorrupt:
      estr = "improper resource content"; break;
   }

   if ( label ) fprintf( stderr, "%s: ", label );
   fprintf( stderr, "%s\n", estr );
   return errcode;
}

/* --- macho data representations ------------------------------------------ */

double	machotime( int time )
{
   // (time_t) 694310400 = 2 Jan 1992 at 00:00:00 GMT
   
   return time == INVALID_LONG ? NaN : (time - 694310400) / 86400.0;
}

int	timemacho( double days )
{
   return isnan(days) ? INVALID_LONG : int(days * 86400 + 694310400.5);
}

int     macho_ra_string( const char *str )
{
   const int precision = (int) RA_DEC_PRECISION;

   int hrs, min;
   float sec;
   if ( sscanf( str, "%d%*c%d%*c%f", &hrs, &min, &sec ) != 3 ) return 0;

   return ( hrs * 15 * 3600 * precision + min * 15 * 60 * precision
	   + int( 0.5 + sec * 15 * precision )
	   );
}
 
char *	macho_ra_string( int ra, char delim, int prcsn )
{
   static char str[64];
   const int precision = (int) RA_DEC_PRECISION;

   int hrs = ra / ( 15 * 3600 * precision );
   ra -= hrs * ( 15 * 3600 * precision );

   int min = ra / ( 15 * 60 * precision );
   ra -= min * ( 15 * 60 * precision );

   float sec = ra / ( 15.0 * precision );

   int width = 3 + prcsn;
   sprintf( str, "%02d%c%02d%c%0*.*f"
	   , hrs, delim, min, delim, width, prcsn, sec );
   return str;
}

int     macho_dec_string( const char *str )
{
   const int precision = (int) RA_DEC_PRECISION;

   int negative;
   while ( isspace(*str) ) ++str;
   if ( negative = (*str == '-') ) ++str;

   int deg, min;
   float sec;
   if ( sscanf( str, "%d%*c%d%*c%f", &deg, &min, &sec ) != 3 ) return 0;

   int dec = ( deg * 3600 * precision + min * 60 * precision
	      + int( 0.5 + sec * precision )
	      );
   return negative ? -dec : dec;
}
 
char *	macho_dec_string( int dec, char delim, int prcsn )
{
   static char str[64];
   const int precision = (int) RA_DEC_PRECISION;

   int negative = (dec < 0);
   if ( negative ) dec = -dec;

   int deg = dec / ( 3600 * precision );
   dec -= deg * ( 3600 * precision );

   int min = dec / ( 60 * precision );
   dec -= min * ( 60 * precision );

   float sec = dec / float(precision);
   
   int width = 3 + prcsn;
   sprintf( str, "%c%02d%c%02d%c%0*.*f"
	   , negative?'-':'+', deg, delim, min, delim, width, prcsn, sec );
   return str;
}

/* --- ChunkToken methods -------------------------------------------------- */

ChunkToken::ChunkToken( SideOfPier SOP, Color CLR, int ROW, int COL )
{
   tok = -1;
   if ( ROW < 0  ||  ROW > 7  ||  COL < 0  ||  COL > 7 ) return;
   tok = (SOP==West ? 0200 : 0) | (CLR==Blue ? 0100 : 0) | (ROW << 3) | COL;
}

ChunkToken::ChunkToken( SideOfPier SOP, short ID )
{
   tok = -1;
   for ( int r = 0; r < 8; ++r )
      for ( int c = 0; c < 8; ++c )
	 if ( ChunkToken(SOP,Red,r,c).id() == ID ) {
	    tok = ChunkToken(SOP,Red,r,c).tok;
	    return;
	 }
	 else if ( ChunkToken(SOP,Blue,r,c).id() == ID ) {
	    tok = ChunkToken(SOP,Blue,r,c).tok;
	    return;
	 }
}

ChunkToken::ChunkToken( SideOfPier SOP, short chip, int ROW, int COL )
{
   tok = -1;
   if ( chip<0 || chip>127 || ROW<0 || ROW>3 || COL<0 || COL>3 ) return;
   short ID = (chip * 16 + 120) % 128;
   tok = ChunkToken( SOP, ID ).tok;
   tok = (tok & 0xFFE4) | (ROW << 3) | COL;
}

short		ChunkToken::chip( void ) const
{
   if ( ! *this ) return -1;
   return ( (8+id()) % 128 ) / 16;
}

short		ChunkToken::id( void ) const
{
   if ( ! *this ) return -1;
   int r = sop()==East ? row() : 7 - row();
   int c = sop()==East ? col() : 7 - col();
   return clr()==Red ? red[r][c] : blue[r][c];
}

ChunkToken	ChunkToken::operator+ ( SideOfPier SOP ) const
{
   return sop() == SOP ? ChunkToken(tok) : -(*this);
}

ChunkToken	ChunkToken::operator+ ( Color CLR ) const
{
   return clr() == CLR ? ChunkToken(tok) : ~(*this);
}

ChunkToken	ChunkToken::incr( int ROW, int COL ) const
{
   return !*this ? ChunkToken() : ChunkToken(sop(),clr(),row()+ROW,col()+COL);
}

const short ChunkToken::red[8][8] = {
   { 12,   13,   14,   15,    0,    4,  120,  124 },
   {  8,    9,   10,   11,    1,    5,  121,  125 },
   { 20,   21,   22,   23,    2,    6,  122,  126 },
   { 16,   17,   18,   19,    3,    7,  123,  127 },
   { 28,   29,   30,   31,   51,   50,   49,   48 },
   { 24,   25,   26,   27,   55,   54,   53,   52 },
   { 36,   37,   38,   39,   43,   42,   41,   40 },
   { 32,   33,   34,   35,   47,   46,   45,   44 }
};
   
const short ChunkToken::blue[8][8] = {
   { 92,   88,  100,   96,  111,  110,  109,  108 },
   { 93,   89,  101,   97,  107,  106,  105,  104 },
   { 94,   90,  102,   98,  119,  118,  117,  116 },
   { 95,   91,  103,   99,  115,  114,  113,  112 },
   { 83,   87,   75,   79,   67,   71,   59,   63 },
   { 82,   86,   74,   78,   66,   70,   58,   62 },
   { 81,   85,   73,   77,   65,   69,   57,   61 },
   { 80,   84,   72,   76,   64,   68,   56,   60 }
};

/* --- unix utilities ------------------------------------------------------ */

double	juliantime( int time )
{
   // (Julian Date) 2440587.5 = (time_t) 0

   return (time / 86400.0) + 2440587.5;
}

int	timejulian( double jd )
{
   return (int) floor( (jd - 2440587.5) * 86400 + 0.5 );
}

int	signal_catch( int sig, sigf_t handler, int flags )
{
   struct sigaction sa;
   memset( &sa, 0, sizeof(sa) );
   *((sigf_t*)&sa.sa_handler) = handler;
   sigemptyset(&sa.sa_mask);
   sa.sa_flags = flags;
   return sigaction(sig,&sa,0);
}

int	signal_unblock( int sig )
{
   sigset_t set;
   sigemptyset(&set);
   if ( sigaddset(&set,sig) ) return -1;
   return sigprocmask(SIG_UNBLOCK,&set,0);
}

int	signal_block( int sig )
{
   sigset_t set;
   sigemptyset(&set);
   if ( sigaddset(&set,sig) ) return -1;
   return sigprocmask(SIG_BLOCK,&set,0);
}

#ifdef __sunos4__

char *	strerror( int i )
{
   extern int sys_nerr;
   extern char *sys_errlist[];
   return ( i >= 0  &&  i < sys_nerr ) ? sys_errlist[i] : NULL;
}

#endif

/* --- lcdir support ------------------------------------------------------- */

//	if buf is a null-pointer the returned string is allocated with new
//	otherwise buf is assumed to have room for MAXPATHLEN+1 bytes

char*
lcdir_star_dir( const char *lcdir, const char *starid, int mode, char *buf )
{
   if ( !lcdir || !starid || !*starid ) return 0;
   int max = strlen( lcdir ) + strlen( starid ) + 8;
   if ( max > MAXPATHLEN ) return 0;
   if ( !buf ) buf = new char [ max+1 ];
   if ( !buf ) return 0;

   const char *hash = starid + strlen(starid) - 2;
   if ( hash < starid ) hash = starid;

   strcpy( buf, lcdir );
   if ( mode == 'w' && mkdir(buf,0777) == -1 && errno != EEXIST ) return 0;
   strcat( buf, "/H" );
   if ( mode == 'w' && mkdir(buf,0777) == -1 && errno != EEXIST ) return 0;
   sprintf( buf + strlen(buf), "/h%s", hash );
   if ( mode == 'w' && mkdir(buf,0777) == -1 && errno != EEXIST ) return 0;
   sprintf( buf + strlen(buf), "/%s", starid );
   if ( mode == 'w' && mkdir(buf,0777) == -1 && errno != EEXIST ) return 0;

   return buf;
}

char*
lcdir_sodlc_file( const char *lcdir, const char *starid, int ext, char *buf )
{
   if ( !lcdir || !starid ) return 0;
   int max = strlen( lcdir ) + strlen( starid ) + 16;
   if ( max > MAXPATHLEN ) return 0;
   if ( !buf ) buf = new char [ max+1 ];
   if ( !buf ) return 0;

   if ( ! lcdir_star_dir( lcdir, starid, 'r', buf ) ) return 0;
   sprintf( buf + strlen(buf), "/sodlc.%c", ext );
   return buf;
}

int lcdir_sodlc_find( const char *lcdir, const char *starid, int source
		      , class SodophotSet &ss )
{
   if ( !lcdir || !starid ) return ss.errcode = SiteErr::errParams;

   int field, index;
   unsigned short tile, seqn;
   char fname[ MAXPATHLEN ];

   if ( sscanf( starid, "%d.%hu.%hu", &field, &tile, &seqn ) != 3 )
      if ( source == 'g' && sscanf( starid, "%d.%d", &field, &index ) == 2 ) {
	 field /= 1000;
	 tile = (unsigned short) (index >> 15);
	 seqn = (unsigned short) (index & 32767);
      }
      else return ss.errcode = SiteErr::errParams;

   if ( lcdir_sodlc_file( lcdir, starid, source, fname ) && !ss.open( fname )
	&& ss.hdr->source == source && ss.hdr->tile == tile
	&& ss.hdr->field == field && ss.star[0].seqn == seqn )
      return 0;

   if ( strchr( "np", source ) ) {
      SodophotSet ss0;
      sprintf( fname, "%s/PhotSS/F_%d/sodset_%d.%d"
	       , lcdir, field, field, tile );
      if ( !ss0.open( fname )
	   && ss0.hdr->source == source && ss0.hdr->tile == tile
	   && ss0.hdr->field == field && (index = ss0.find(seqn)) != -1 )
	 return excerpt( ss, ss0, 1, &index, 0, 0 );
   }

   return ss.errcode = SiteErr::errMissing;
}

/* --- Initialize IEEE Infinity and NaN ------------------------------------ */

double Infinity, NaN;

void initialize_site_globals( void )
{
   if ( isnan(NaN) ) return;	// initialize just once

#if defined __sunos4__ || defined __solaris__ || defined _IEEE_FP

   double Zero = 0;
   Infinity = 1 / Zero;
   NaN = Infinity - Infinity;

#elif defined __GNUG__

   Infinity = 1.8e308;
   NaN = 1.8e308 - 1.8e308;

#else

   float inf, nan;
   *((int*) &inf) = 0x7F800000;
   *((int*) &nan) = 0x7FC00000;
   Infinity = inf;
   NaN = nan;

#endif

   if ( ! isnan(NaN) ) {
      fprintf( stderr, "initialize_site_globals: "
	       "failed to generate IEEE NaN, aborting\n" );
      exit( 1 );
   }
}

class InitSiteGlobals {
public:
   InitSiteGlobals( void ) { initialize_site_globals(); }
};
static InitSiteGlobals site;
