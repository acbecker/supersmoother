/* written by John Doug Reynolds, March 1997 */

#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "fits_header.h"

#define FITS_HEADER_BYTES 2880
#define FITS_HEADER_CARDS 36

static	char* h_key_END = "END";
static	char* fix_key( char* str );
static  char* floating( char* type, void* value );
static  int is_blank( HEADER hdr, int card );
static	void rubout( HEADER hdr, int card1, int card2 );
static	int add_pages( HEADER hdr, int new_end, int erase );

/* convert error codes into friendly messages */
char* h_strerror( int error_code )
{
   switch ( error_code ) {
   case h_no_err:
      return "fits_header: no error";
   case h_logic_err:
      return "fits_header: request is illogical or illegal";
   case h_format_err:
      return "fits_header: header format is fucked";
   case h_mem_err:
      return "fits_header: memory allocation error";
   case h_file_err:
      return "fits_header: file access error or unexpected eof";
   default:
      return "fits_header: error error (fubar)";
   }
}

/* create a brand new header, properly initialized with just an END card */
HEADER new_header()
{
   HEADER hdr = (HEADER) malloc( sizeof(struct header_control) );
    
   if ( hdr ) {
      hdr->end = 0;
      hdr->H = (CARD*) malloc( FITS_HEADER_BYTES );
      if ( !hdr->H ) {
	 free( hdr );
	 return NULL;
      }
      rubout( hdr, 0, FITS_HEADER_CARDS - 1 );
      h_write_card( hdr, hdr->end, h_key_END, NULL, NULL, NULL );
   }
   return hdr;
}

/* deallocate the space used by the header */
void free_header( HEADER hdr )
{
   if ( hdr ) {
      free( hdr->H );
      free( hdr );
   }
}

/* return total size of header in bytes, or -1 on error */
off_t header_offset( HEADER hdr )
{
   if ( !hdr ) return (off_t) -1;
   return (off_t) ( (1 + hdr->end / FITS_HEADER_CARDS) * FITS_HEADER_BYTES );
}

/* add the specified number of blank cards to the end of the header */
int enlarge_header( HEADER hdr, int cards )
{
   int	old_end, error;

   if ( !hdr ) return h_logic_err;

   old_end = hdr->end;
   if ( error = add_pages( hdr, old_end + cards, 1 ) ) return error;
   rubout( hdr, old_end, hdr->end - 1 );
   return h_write_card( hdr, hdr->end, h_key_END, NULL, NULL, NULL );
}

/* create a new header containing a copy of the given header */
HEADER copy_header( HEADER srce_hdr )
{
   HEADER	dest_hdr;
   off_t	hdr_bytes;

   if ( (hdr_bytes = header_offset( srce_hdr )) == -1 )
      return NULL;
   if ( !(dest_hdr = (HEADER) malloc( sizeof(struct header_control) )) )
      return NULL;
   if ( !(dest_hdr->H = (CARD*) malloc( hdr_bytes )) ) {
      free( dest_hdr );
      return NULL;
   }
   memcpy( dest_hdr->H, srce_hdr->H, hdr_bytes );
   dest_hdr->end = srce_hdr->end;
   return dest_hdr;
}

/*
o  write header to file only if END card can be found
o  if END card is found but not where expected, treachery is suspected
      and a warning is sent to stderr; hdr->end is reset to the card where
      END was actually found and remaining header pages are not written;
      this behavior protects the data from being confused with the header
o  this function does not free the header after writing it
*/   
static int u_write_header( FILE* fp, int fd, HEADER hdr )
	/* if fp is NULL then fd must be non-negative */
{
   off_t	hdr_bytes;
   off_t	bytes_sent;
   int		card = 0;
    
   if ( !hdr ) return h_logic_err;
   if ( !h_find_card(hdr,h_key_END,&card) ) return h_format_err;
    
   if ( card != hdr->end ) {
      fprintf( stderr, "fits_header: warning, END card out of place\n" );
      hdr->end = card;
   }

   hdr_bytes = header_offset( hdr );
    
   if ( fp )
      bytes_sent = fwrite( hdr->H, 1, hdr_bytes, fp );
   else
      bytes_sent = write( fd, hdr->H, hdr_bytes );

   return ( bytes_sent == hdr_bytes ) ? h_no_err : h_file_err;
}
int fwrite_header( FILE* fi, HEADER hdr )
{
   return u_write_header( fi, -1, hdr );
}
int write_header( int fd, HEADER hdr )
{
   return u_write_header( NULL, fd, hdr );
}

/*
o  create a new header by reading page at a time from the given file
	until a page is read which contains the END card
o  returns NULL if end of file is encountered before the END card
o  if header_size is non-null, *header_size will receive the length
	of the header in bytes
o  if err is non-null, *err will receive the error status code
*/
static HEADER u_read_header( FILE* fp, int fd, size_t* header_size, int* err )
	/* if fp is NULL then fd must be non-negative */
{
   HEADER	hdr	= new_header();
   int		error	= h_no_err;
   int		card	= 0;

   if ( !hdr ) {
      if ( err ) *err = h_mem_err;
      return NULL;
   }
   else hdr->end = FITS_HEADER_CARDS - 1;

   while ( !error ) {
      off_t bytes_read;
	
      if ( fp )
	 bytes_read = fread( hdr->H[card], 1, FITS_HEADER_BYTES, fp );
      else
	 bytes_read = read( fd, hdr->H[card], FITS_HEADER_BYTES );
      
      if ( bytes_read != FITS_HEADER_BYTES ) {
	 if ( bytes_read == -1 ) perror( NULL );
	 error = h_file_err;
	 break;
      }

      if ( error ) break;
	
      if ( h_find_card( hdr, h_key_END, &card ) ) {
	 hdr->end = card;
	 if ( header_size ) *header_size = header_offset( hdr );
	 if ( err ) *err = h_no_err;
	 return hdr;
      }
      else {
	 card += FITS_HEADER_CARDS;
	 error = add_pages( hdr, hdr->end + FITS_HEADER_CARDS, 0 );
      }
   }

   free_header( hdr );
   if ( header_size ) *header_size = 0;
   if ( err ) *err = error;
   return NULL;
}
HEADER fread_header( FILE* fi, size_t* header_size, int* err )
{
   return u_read_header( fi, -1, header_size, err );
}
HEADER read_header( int fd, size_t* header_size, int* err )
{
   return u_read_header( NULL, fd, header_size, err );
}

char* h_get_card( HEADER hdr, int card )
{
    if ( !hdr  ||  card < 0  ||  card > hdr->end ) return NULL;

    return hdr->H[card];
}

char* h_find_card( HEADER hdr, char* key, int* card_num )
{
    int	    card   = (card_num && *card_num >= 0) ? *card_num : 0;
    char*   key8   = fix_key( key );

    if ( !key  ||  !hdr  ||  card > hdr->end ) return NULL;

    while ( strncmp( key8, hdr->H[card], 8 ) != 0 )
       if ( ++card > hdr->end ) return NULL;

    if ( card_num ) *card_num = card;
    return hdr->H[card] + 10;
}

char* h_card_search( HEADER hdr, char* str, int* card_num )
{
    int   length = strlen( str );
    int   card   = (card_num && *card_num >= 0) ? *card_num : 0;
    char* pc;

    if ( !hdr  ||  length == 0 ) return NULL;

    for (; card <= hdr->end; ++card ) {
	int count = 80 - length + 1;
	
	for ( pc = hdr->H[card]; count--; ++pc )
	    if ( strncmp(pc,str,length) == 0 ) {
		if ( card_num ) *card_num = card;
		return pc + length;
	    }
    }
    return NULL;
}

int h_read_card( HEADER hdr, char* key, char* type, void* value )
{
   char str[80], *ptr;
   char* pc = h_find_card( hdr, key, NULL );
   char* start = pc - 10;

   /* special treatment for type "B"; *value = (card exists) */
   if ( *type == 'B' ) {
      *( (int*) value ) = pc ? 1 : 0;
      return h_no_err;
   }

   if (!pc) return h_format_err;
   if ( !type  ||  !value ) return h_logic_err;
    
   switch ( *type ) {
   case 'b':
      *( (int*) value ) = (pc[19] == 'T');
      break;
   case 's':
      ptr = (char*) value;
      for ( ++pc; *pc != '\047' && pc-start < 79; *ptr++ = *pc++ );
      while ( ptr > (char*)value  &&  *(ptr-1) == ' ' ) --ptr;
      *ptr = '\0';
      break;
   default:
      str[0] = '%'; str[1] = '\0';
      strcat( str, type );
      if ( !sscanf( pc, str, value ) ) return 2;
   }

   return h_no_err;
}

/* data size modifier 'h' applies only to type 'd'
                      'L' applies only to type 'f'
		      'l' applies only to types 'd' and 'f'
   if key is NULL, the keyword field is not overwritten
   if value is NULL, the value field is blank with no '='
   if comment is NULL, the comment field is blank with no '/'
*/
int h_write_card( HEADER hdr, int card, char* key,
		 char* type, void* value, char* comment )
{
   char *pc, *start, *ptr;

   if ( !(start = h_get_card(hdr,card)) ) return h_logic_err;

   if ( key ) strcpy( start, fix_key(key) );

   /* write in the equals sign if there is a value */
   pc = start + 8;
   *pc++ = ( value ) ? '=' : ' ';
   *pc++ = ' ';

   /* write the formatted value */

   if ( type  &&  value ) {
      if ( strcmp( type, "b" ) == 0 ) {
	 sprintf( pc, "%20s", *((int*)value) ? "T" : "F" );
	 pc += 20;
      }
      else if ( strcmp( type, "hd" ) == 0 ) {
	 sprintf( pc, "%20hd", *((short*)value) );
	 pc += 20;
      }
      else if ( strcmp( type, "d" ) == 0 ) {
	 sprintf( pc, "%20d", *((int*)value) );
	 pc += 20;
      }
      else if ( strcmp( type, "ld" ) == 0 ) {
	 sprintf( pc, "%20ld", *((long*)value) );
	 pc += 20;
      }
      else if ( strchr( type, 'f' )  ||  strchr( type, 'e' ) ) {
	 sprintf( pc, "%20s", floating(type,value) );
	 pc += 20;
      }
      else if ( strcmp( type, "s" ) == 0 ) {
	 *pc++ = '\047';
	 for ( ptr = (char*)value; *ptr && pc-start < 79; *pc++ = *ptr++ );
	 while ( pc-start < 19 ) *pc++ = ' ';
	 *pc++ = '\047';
      }
      else if ( strcmp( type, " " ) == 0 ) {
	 /* hanging comment; nothing is required here */
      }
      else
	 return h_logic_err;
   }
    
   if ( comment  &&  pc - start < 70 ) {
      if ( type ) {
	 while ( pc - start < 40 ) *pc++ = ' ';
	 *pc++ = '/'; *pc++ = ' ';
      }
      while ( *comment  &&  pc - start < 80 ) *pc++ = *comment++;
   }
   while ( pc - start < 80 ) *pc++ = ' ';
    
   if ( comment  &&  *comment ) {
      if ( type ) type = " ";
      if ( is_blank( hdr, card+1 ) )
	 return h_write_card( hdr, card+1, "COMMENT", type, NULL, comment );
      else if ( card+1 == hdr->end )
	 return h_append_card( hdr, "COMMENT", type, NULL, comment );
   }
   return h_no_err;
}

/*
   Write the card at the end of the header.  If a block of blank lines is
   found just before the END card, the first blank line in the block will
   be overwritten.  Otherwise, the card is appended.
*/
int h_append_card( HEADER hdr, char* key,
		   char* type, void* value, char* comment )
{
   int card, error;
    
   if ( !hdr  ||  !key ) return h_logic_err;
    
   /* starting from END card search back to the first non-blank card */
   for ( card = hdr->end - 1; is_blank( hdr, card ); --card );

   if ( card < hdr->end - 1 )
      ++card;
   else {
      if ( error = enlarge_header( hdr, 1 ) ) return error;
      card = hdr->end - 1;
   }

   return h_write_card( hdr, card, key, type, value, comment );
}

int h_change_card( HEADER hdr, char* key,
		   char* type, void* value, char* comment )
{
   int card = 0;
    
   if ( !hdr  ||  !key ) return h_logic_err;

   if ( h_find_card( hdr, key, &card ) )
      return h_write_card( hdr, card, key, type, value, comment );
   else
      return h_append_card( hdr, key, type, value, comment );
}

/* ------------------------------------------------------------------------- */
/* ------------------------  proprietary functions  ------------------------ */

/*
    add_pages sets hdr->end to the specified value and creates enough
    new pages to make hdr->end valid.  initializing the new space with
    space characters is optional.  hdr->end may not be moved backward.
    does not erase the old END card or create a new one.
*/
static int add_pages( HEADER hdr, int new_end, int erase )
{
   size_t old_pages, new_pages;

   if ( !hdr  ||  hdr->end >= new_end ) return h_logic_err;

   old_pages = 1 + hdr->end / FITS_HEADER_CARDS;
   new_pages = 1 + new_end / FITS_HEADER_CARDS;

   if ( new_pages > old_pages ) {
      void* temp = realloc( hdr->H, new_pages * FITS_HEADER_BYTES );
      if ( !temp ) return h_mem_err;
      hdr->H = (CARD*) temp;
   }

   if ( erase )
      rubout(hdr,old_pages*FITS_HEADER_CARDS,new_pages*FITS_HEADER_CARDS-1);

   hdr->end = new_end;
   return h_no_err;
}

/*
   rubout should be used with care.  rubout writes space characters over
   all the cards from card1 to card2, inclusive.
   rubout does not verify that card2 < hdr->end.
*/
static void rubout( HEADER hdr, int card1, int card2 )
{
   char*	pc = hdr->H[card1];
   int		count	= 80 * ( card2 - card1 + 1 );

   while ( count-- ) *pc++ = ' ';
}

/*
   fix_key returns a null-terminated eight character version of str
   in a static nine character array
*/
static char* fix_key( char* str )
{
   static char key[9];
   int n;
   if ( !str ) str = "";
   for ( n = 0; *str  &&  n < 8; ++n ) key[n] = *str++;
   while ( n < 8 ) key[n++] = ' ';
   key[8] = '\0';
   return key;
}

/*
   is_blank returns true if the 80 character array pc is blank
*/
static int is_blank( HEADER hdr, int card )
{
   int count;
   char *pc = h_get_card( hdr, card );
   if ( !pc ) return 0;
   for ( count = 80; count--; ++pc ) if ( *pc && *pc != ' ' ) break;
   return (count < 0);
}

/*
   floating formats floating point numbers with the decimal point and
   at least one following digit; within this constraint, trailing zeros
   are removed.  This approach has aesthetic appeal and should be
   acceptable for the screwy MACHO version of IRAF.
*/
static char* floating( char* type, void* value )
{
   static char	str[80];
   char		len, chr, *fmt, *srce, *dest;

   if ( !type  ||  !*type  ||  strlen(type) > (size_t)2 ) return NULL;
   else if ( strlen(type) == 1 ) { len = 'h'; chr = *type; }
   else { len = type[0]; chr = type[1]; }

   if ( chr == 'f' ) fmt = "%#.*G";
   else if ( chr == 'e' ) fmt = "%#.*E";
   else return NULL;
   
   if ( len == 'h' )
      sprintf( str, fmt, 7, *((float*)value) );
   else if ( len == 'l' )
      sprintf( str, fmt, 13, *((double*)value) );
   else if ( len == 'L' )
      sprintf( str, fmt, 13, *((long double*)value) );
   else
      return NULL;

   if ( !(srce = strchr( str, '.' )) ) return str;
   
   for ( dest = srce += 2; isdigit(*srce); ++srce )
      if ( *srce != '0' ) dest = srce + 1;
   while ( *dest++ = *srce++ );
   
   return str;
}
