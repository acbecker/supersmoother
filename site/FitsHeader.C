/* written by John Doug Reynolds, September 1996 */

#include "libsite.H"
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

static const int pagelength = 36;
static const int pagesize = 2880;
static const char end_key[]	= "END     ";
static const char history_key[] = "HISTORY ";
static const char comment_key[] = "COMMENT ";
static const char blank_card[] =
	"                                        "
	"                                        ";

int	FitsHeader::grow( int count )
{
   errcode = 0;

   if ( count < 0 ) return errcode = errParams;
   if ( count == 0 ) return 0;

   typedef char Card[80];
   Card *temp = new Card [ (pages + count) * pagelength ];
   if ( !temp ) return errcode = errMemory;

   if ( card ) {
      memcpy( temp[0], card[0], pages * pagesize );
      delete [] card;
   }
   memset( temp[pages * pagelength], ' ', count * pagesize );

   card = temp;
   pages += count;
   return 0;
}

FitsHeader::FitsHeader( int count )
{
   memset( buffer, 0, sizeof(buffer) );
   cursor = end = -1;
   pages = 0;
   card = 0;
   grow( count );
}

FitsHeader::~FitsHeader( void )
{
   if ( card ) delete [] card;
}

FitsHeader::operator char*() const
{
   return cursor >= 0 && cursor <= end ? card[cursor] : 0;
}

char*	FitsHeader::operator[] ( int index ) const
{
   return index >= 0 && index <= end ? card[index] : 0;
}

int	FitsHeader::length( void ) const
{
   return end + 1;
}

int	FitsHeader::size( void ) const
{
   return ((end + pagelength) / pagelength) * pagesize;
}

char	FitsHeader::classify( int index ) const
{
   if ( index == -1 ) index = cursor;
   if ( index < 0 || index > end ) return 0;
   char *str = card[index];
   if ( strncmp( str, comment_key, 8 ) == 0 ) {
      int i;
      for ( i = 8; i < 80 && str[i] == ' '; ++i );
      return (i < 80 && str[i] == '/') ? '/' : 'C';
   }
   if ( strncmp( str, end_key, 8 ) == 0 ) return 'E';
   if ( strncmp( str, history_key, 8 ) == 0 ) return 'H';
   if ( strncmp( str, blank_card, 80 ) == 0 ) return 'B';
   return '*';
}

FitsHeader& FitsHeader::operator= ( const FitsHeader &hdr )
{
   errcode = 0;
   end = -1;

   if ( hdr.length() > 0 ) {
      int count = hdr.size() / pagesize;
      if ( count > pages && grow(count-pages) ) return *this;
      memcpy( card[0], hdr[0], hdr.size() );
      end = hdr.length() - 1;
   }

   return *this;
}

int	FitsHeader::extend( int count )
{
   errcode = 0;

   // allocate any extra pages required and initialize if necessary

   int total = 1 + (abs(count) + (end < 0 ? 0 : end)) / pagelength;
   if ( total > pages && grow(total-pages) ) return errcode;

   if ( end < 0 ) {
      memcpy( card[0], end_key, 8 );
      if ( count ) cursor = 0;
      end = 0;
   }

   if ( count == 0 ) return 0;

   // move following cards forward and erase the vacated space

   if ( count < 0 )
      count = -count;
   else
      cursor = end;

   if ( cursor < 0 || cursor > end ) return errcode = errParams;

#ifdef __sunos4__
   bcopy( card[cursor], card[cursor+count], 80 * (end-cursor+1) );
#else
   memmove( card[cursor+count], card[cursor], 80 * (end-cursor+1) );
#endif
   memset( card[cursor], ' ', 80 * count );
   end += count;

   return 0;
}

int	FitsHeader::print( int index, FILE *fi )
{
   errcode = 0;
   if ( index == -1 ) index = cursor;
   if ( index < 0 || index > end ) return errcode = errParams;
   fprintf( fi, "%.80s\n", card[index] );
   return 0;
}

int	FitsHeader::write( int fd )
{
   errcode = 0;
   int nbytes = size();
   if ( classify(end) != 'E' ) return errcode = errCorrupt;
   if ( ::write( fd, card[0], nbytes ) != nbytes ) return errcode = -errno;
   return 0;
}

int	FitsHeader::read( int fd )
{
   int i = 0;
   errcode = 0;
   end = -1;

   for (;;) {
      if ( i/pagelength >= pages && grow(1) )
	 return errcode;
      if ( ::read( fd, card[i], pagesize ) != pagesize )
	 return errcode = -errno;
      for ( int m = i + pagelength; i < m; ++i )
	 if ( strncmp( card[i], end_key, 8 ) == 0 ) {
	    end = i;
	    return 0;
	 }
   }
}

FitsHeader&	FitsHeader::key( const char *str, int seek )
{
   if ( !strchr("^$><",seek) || !str ) {
      cursor = -1;
      errcode = errParams;
      return *this;
   }

   int increment = ( seek == '^' || seek == '>' ) ? 1 : -1;
   cursor = (seek != '^') ? (seek != '$') ? cursor + increment : end : 0;
   errcode = 0;

   char key8[9];
   sprintf( key8, "%-8.8s", str );
   for ( char *p; p = *this; cursor += increment )
      if ( strncmp( p, key8, 8 ) == 0 ) return *this;

   cursor = -1;
   errcode = errMissing;
   return *this;
}

FitsHeader&	FitsHeader::search( const char *str, int seek )
{
   int length = str ? strlen( str ) : 0;
   if ( !strchr("^$><",seek) || length <= 0 || length > 70 ) {
      cursor = -1;
      errcode = errParams;
      return *this;
   }

   int increment = ( seek == '^' || seek == '>' ) ? 1 : -1;
   cursor = (seek != '^') ? (seek != '$') ? cursor + increment : end : 0;
   errcode = 0;

   int count, tries = 70 - length + 1;
   for ( char *p = 0; p = *this; cursor += increment )
      for ( count = tries, p += 10; count--; ++p )
	 if ( strncmp(p,str,length) == 0 ) return *this;

   cursor = -1;
   errcode = errMissing;
   return *this;
}

// width() returns the number of characters in the value portion of a
// FITS card.  The value portion is terminated by the first '/' that
// isn't within a FITS-style string.

static int width( const char *str )
{
   int i = 10;
   while ( i < 80 && str[i] == ' ' ) ++i;
   if ( i < 80 && str[i] == '\'' )
      while ( ++i < 80 && str[i] != '\'' );
   while ( i < 80 && str[i] != '/' ) ++i;
   return i - 10;
}

// get() transcribes the value portion of the current card into buffer.
// If there is no current card the previous error is not disturbed if
// non-zero.  Returns buffer or zero on error.

char*	FitsHeader::get( void )
{
   buffer[0] = '\0';
   char *str = *this;
   if ( ! str ) {
      if ( ! errcode ) errcode = errParams;
      return 0;
   }
   if ( str[8] != '=' ) {
      errcode = errFormat;
      return 0;
   }
   int n = width(str);
   memcpy( buffer, str+10, n );
   buffer[n] = '\0';
   errcode = 0;
   return buffer;
}

int	FitsHeader::operator>> ( char &value )
{
   char v[2], *str = get();
   if ( str && sscanf(str,"%1s",v) != 1 ) errcode = errFormat;
   if ( !errcode ) value = v[0];
   return errcode;
}

int	FitsHeader::operator>> ( short &value )
{
   char *str = get();
   if ( str && sscanf(str,"%hi",&value) != 1 ) errcode = errFormat;
   return errcode;
}

int	FitsHeader::operator>> ( int &value )
{
   char *str = get();
   if ( str && sscanf(str,"%i",&value) != 1 ) errcode = errFormat;
   return errcode;
}

int	FitsHeader::operator>> ( long &value )
{
   char *str = get();
   if ( str && sscanf(str,"%li",&value) != 1 ) errcode = errFormat;
   return errcode;
}

int	FitsHeader::operator>> ( float &value )
{
   char *str = get();
   if ( str && sscanf(str,"%f",&value) != 1 ) errcode = errFormat;
   return errcode;
}

int	FitsHeader::operator>> ( double &value )
{
   char *str = get();
   if ( str && sscanf(str,"%lf",&value) != 1 ) errcode = errFormat;
   return errcode;
}

int	FitsHeader::operator>> ( char *value )
{
   char *src = get(), *dst = value;
   if ( errcode ) return errcode;
   if ( *src++ != '\'' ) return errcode = errFormat;
   while ( *src && *src != '\'' ) *dst++ = *src++;
   while ( dst > value && *(dst-1) == ' ' ) --dst;
   *dst = '\0';
   return 0;
}

FitsHeader&	FitsHeader::operator<< ( char value )
{
   errcode = 0;
   sprintf( buffer, "%20c", value );
   return *this;
}

FitsHeader&	FitsHeader::operator<< ( short value )
{
   errcode = 0;
   sprintf( buffer, "%20d", value );
   return *this;
}

FitsHeader&	FitsHeader::operator<< ( int value )
{
   errcode = 0;
   sprintf( buffer, "%20d", value );
   return *this;
}

FitsHeader&	FitsHeader::operator<< ( long value )
{
   errcode = 0;
   sprintf( buffer, "%20d", value );
   return *this;
}

FitsHeader&	FitsHeader::operator<< ( float value )
{
   errcode = 0;
   char temp[32], *src, *dst;
   sprintf( temp, "%#.6G", value );
   for ( dst = src = strchr(temp,'.') + 2; isdigit(*src); ++src )
      if ( *src != '0' ) dst = src + 1;
   while ( *dst++ = *src++ );
   sprintf( buffer, "%20s", temp );
   return *this;
}

FitsHeader&	FitsHeader::operator<< ( double value )
{
   errcode = 0;
   char temp[32], *src, *dst;
   sprintf( temp, "%#.15G", value );
   for ( dst = src = strchr(temp,'.') + 2; isdigit(*src); ++src )
      if ( *src != '0' ) dst = src + 1;
   while ( *dst++ = *src++ );
   sprintf( buffer, "%20s", temp );
   return *this;
}

FitsHeader&	FitsHeader::format( double value,int conversion,int precision )
{
   errcode = 0;
   char fstr[18];
   sprintf( fstr, "%%#20.%d%c", precision, conversion );
#ifdef __sunos4__
   if ( ! sprintf( buffer, fstr, value ) ) errcode = errParams;
#else
   if ( sprintf( buffer, fstr, value ) < 0 ) errcode = errParams;
#endif
   return *this;
}

FitsHeader&	FitsHeader::operator<< ( const char *value )
{
   errcode = 0;
   sprintf( buffer, "'%-8.68s'", value );
   return *this;
}

int	FitsHeader::insert( const char *key, const char *rem, const char *val )
{
   if ( val == (char*) -1 ) val = buffer;
   int vlen = val ? strlen(val) : 0;
   int rlen = rem ? strlen(rem) : 0;

   if ( ! key || vlen > 70 ) return errcode = errParams;

   char key8[9];
   sprintf( key8, "%-8.8s", key );
   key = key8;

   // insert enough extra space for the card and any hanging comments

   int rlen1, lines;

   if ( val ) {
      rlen1 = 68 - (vlen+1 > 30 ? vlen+1 : 30);
      if ( rlen1 < 0 ) rlen1 = 0;
      if ( rlen < rlen1 ) rlen1 = rlen;
      lines = 1 + (rlen + 37 - rlen1) / 38;
   }
   else {
      rlen1 = rlen < 70 ? rlen : 70;
      lines = 1 + (rlen + 69 - rlen1) / 70;
   }

   for ( int i = cursor; lines && classify(i++) == 'B'; --lines );
   if ( extend(-lines) ) return errcode;

   // construct the first card

   char *str = *this;
   memcpy( str, key, 8 );
   if ( val ) str[8] = '=';
   str += 10;
   if ( val ) {
      memcpy( str, val, vlen );
      str += vlen;
   }
   if ( val && rlen1 ) {
      if ( vlen < 29 ) str += 29 - vlen;
      memcpy( str, " / ", 3 );
      str += 3;
   }
   if ( rlen1 ) {
      memcpy( str, rem, rlen1 );
      rlen -= rlen1;
      rem += rlen1;
   }
   ++cursor;

   // construct comment overflow cards

   if ( val ) key = comment_key;
   int width = val ? 38 : 70;

   while ( rlen > 0 ) {
      memcpy( str = *this, key, 8 );
      str += 10;
      if ( val ) {
	 str[30] = '/';
	 str += 32;
      }
      rlen1 = rlen < width ? rlen : width;
      memcpy( str, rem, rlen1 );
      rlen -= rlen1;
      rem += rlen1;
      ++cursor;
   }

   return 0;
}

int	FitsHeader::append( const char *key, const char *rem, const char *val )
{
   for ( cursor = end; cursor > 0 && classify(cursor-1) == 'B'; --cursor );
   return insert(key,rem,val);
}

int	FitsHeader::update( const char *key, const char *rem, const char *val )
{
   char temp[73];

   // if the key does not already exist, just append this card

   if ( ! key ) return errcode = errParams;
   if ( this->key(key).errcode ) return append(key,rem,val);

   // if there is a new remark, erase old hanging comments

   if ( rem )
      for ( int i = cursor + 1; classify(i) == '/'; ++i )
	 memset( card[i], ' ', 80 );

   // otherwise create a new remark from the original

   else {
      char *str = card[cursor];
      int j = 79, i = (str[8] == '=') ? width(str)+11 : 8;
      while ( i < 80 && str[i] == ' ' ) ++i;
      while ( j >= i && str[j] == ' ' ) --j;
      int n = i <= j ? j - i + 1 : 0;
      memcpy( temp, str+i, n );
      temp[n] = '\0';
      rem = temp;
   }

   // erase the old card and replace it with the new one

   memset( card[cursor], ' ', 80 );
   insert(key,rem,val);

   if ( rem == temp )
      while ( classify(cursor) == '/' ) ++cursor;

   return errcode;
}
