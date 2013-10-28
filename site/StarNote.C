/* written by John Doug Reynolds, January 1998 */

#include <string.h>
#include "Macho.h"
#include "libsite.H"

/* --- comments ------------------------------------------------------------

  A StarNote is basically an arbitrary array of bytes, with a thirty
  character label that identifies the data format.  If the note is not
  formatted as a printable, null-terminated ascii string the label
  should have # as its first character. */

/* --- StarNote methods ---------------------------------------------------- */

StarNote::StarNote( void )
{
   buflen = 0;
   field = 0;
   tile = INVALID_USHORT;
   seqn = INVALID_USHORT;
   length = 0;
   buffer = 0;
   label[0] = '\0';
   memset( user, 0, sizeof(user) );
   date = 0;
}

StarNote& StarNote::operator=( const StarNote &note )
{
   if ( set_length(note.length) ) return *this;
   field = note.field;
   tile = note.tile;
   seqn = note.seqn;
   date = note.date;
   strcpy( label, note.label );
   memcpy( user, note.user, 8 );
   if ( note.length ) memcpy( buffer, note.buffer, note.length );
   return *this;
}

StarNote::~StarNote( void )
{
   if ( buffer ) delete [] buffer;
}

int	StarNote::set_length( int len )
{
   if ( len < 0 || len > 65536 ) return 1;

   if ( len <= buflen ) {
      length = (unsigned short) len;
      return 0;
   }

   // extend to nearest multiple of 128
   int nbytes = len + 128 - (len & 127);
   char *newbuf = new char [nbytes];
   if ( !newbuf ) return 1;

   if ( buffer ) {
      memcpy( newbuf, buffer, length );
      delete [] buffer;
   }

   buflen = nbytes;
   buffer = newbuf;
   length = len;
   return 0;
}

int	StarNote::set_label( const char *lbl )
{
   if ( !lbl ) {
      label[0] = '\0';
      return 0;
   }

   if ( strlen(lbl) + 1 > sizeof(label) ) return 1;
   strcpy( label, lbl );
   return 0;
}

int	StarNote::set_user( void )
{
   user[0] = '\0';
   char str[L_cuserid];
   if ( !cuserid(str) ) return 1;
   strncpy( user, str, sizeof(user) - 1 );
   return 0;
}

int	StarNote::print( FILE *fi ) const
{
   fprintf( fi, " %s  %.4f  %d.%d.%d  %s\n"
	    , user, machotime(date), field, tile, seqn, label );
   fprintf( fi,
	    "---------------------------------------"
	    "---------------------------------------\n" );
   if ( label[0] != '#' && length > 1 ) fprintf( fi, "   %s\n", buffer );
   fprintf( fi, "\n" );
   return 0;
}

int	StarNote::scan( FILE *fi )
{
   float fdate;
   char scratch[1024];
   if ( 6 != fscanf( fi, "%s%f%hd.%hu.%hu%s"
		     , user, &fdate, &field, &tile, &seqn, label ) )
      return 1;

   set_length(1);
   buffer[0] = '\0';
   date = (int) timemacho(fdate);

   // discard rest of this line and next line
   if ( !fgets(scratch,1024,fi) || !fgets(scratch,1024,fi) ) return 1;

   if ( label[0] != '#' ) {
      if ( !fgets(scratch,1024,fi) ) return 1;
      int len = strlen(scratch);
      if ( len > 3 ) {
	 if ( scratch[len-1] != '\n' ) return 1;
	 set_length(len-3);
	 memcpy( buffer, scratch+3, length-1 );
	 buffer[length-1] = '\0';
      }
   }
   return 0;
}

/* --- StarNote construction utilities ------------------------------------- */

const char* category_source_menu[] = {
   "unanalyzed", "consensus", "phased",
   "sodstat varcomb", "sodstat microcomb", "multiple alert",
   0
};

const char* category_name_menu[] = {
   "generic variable", "egregious variable", "long-period variable",
   "eruptive variable", "microlensing", "pseudo-microlensing",
   "supernova echo", "eclipsing binary", "Cepheid", "RR Lyrae",
   0
};

const char* startag_alert_menu[] = {
   "AlwaysIgnore", "StandardCuts", "MinimalCuts",
   0
};

const char* select_from_menu( const char *label, const char *list[], FILE *fi )
{
   int items = 0;
   while ( list[items] ) ++items;
   if ( items == 0 ) return 0;

   fprintf( fi, "%s:\n", label );

   int rows = int( 0.5 * items + 0.5 );
   for ( int i = 0; i < rows; ++i ) {
      fprintf( fi, "%4d] %-36s", i+1, list[i] );
      if ( i+rows < items )
	 fprintf( fi, "%4d] %s", i+rows+1, list[i+rows] );
      fprintf( fi, "\n" );
   }

   while ( 1 ) {
      int i;
      char buf[80];
      fprintf( fi, "select 'q' or choice by number [1]: " );
      if ( !fgets(buf,80,stdin) ) { fprintf(fi,"\n"); return 0; }
      if ( strcmp(buf,"q\n") == 0 ) return 0;
      if ( buf[0] == '\n' ) return list[0];
      if ( sscanf(buf,"%d",&i) == 1 && i > 0 && i <= items ) return list[i-1];
      fprintf( fi, "INVALID SELECTION\n" );
   }
}

int construct_category_note( StarNote &note, short field, ushort tile
			     , ushort seqn, const char *src, const char *cat )
{
   if ( !src )
      if ( !(src=select_from_menu("Type of analysis",category_source_menu)) )
	 return 1;

   if ( !cat )
      if ( !(cat=select_from_menu("Type of star",category_name_menu)) )
	 return 1;

   note.set_label( "Category" );
   note.field = field;
   note.tile = tile;
   note.seqn = seqn;
   note.set_user();
   note.date = (int) time(0);

   if ( note.set_length( strlen(src) + strlen(cat) + 4 ) ) return 1;
   sprintf( note.buffer, "%s / %s", src, cat );
   return 0;
}

int construct_startag_note( StarNote &note, short field, ushort tile
			    , ushort seqn, const char *alert, int cameos )
{
   if ( !alert )
      if ( !(alert=select_from_menu("Alert disposition",startag_alert_menu)) )
	 return 1;
   if ( strcmp( alert, "StandardCuts" ) == 0 ) alert = "";
   if ( *alert == '\0' ) alert = 0;

   if ( cameos == -1 )
      while ( 1 ) {
	 char buf[80];
	 printf( "\nAutomatically collect cameos? [n]: " );
	 if ( !fgets(buf,80,stdin) ) { printf("\n"); return 1; }
	 if ( buf[0] == '\n' ) { cameos = 0; break; }
	 if ( strcmp(buf,"n\n") == 0 ) { cameos = 0; break; }
	 if ( strcmp(buf,"y\n") == 0 ) { cameos = 1; break; }
	 printf( "INVALID SELECTION\n" );
      }

   note.set_label( "StarTag" );
   note.field = field;
   note.tile = tile;
   note.seqn = seqn;
   note.set_user();
   note.date = (int) time(0);

   if ( note.set_length(1) ) return 1;
   note.buffer[0] = '\0';

   if ( alert ) {
      if ( note.set_length( note.length + strlen(alert) ) ) return 1;
      strcat( note.buffer, alert );
   }
   if ( alert && cameos ) {
      if ( note.set_length( note.length + 1 ) ) return 1;
      strcat( note.buffer, " " );
   }
   if ( cameos ) {
      const char *str = "CollectCameos";
      if ( note.set_length( note.length + strlen(str) ) ) return 1;
      strcat( note.buffer, str );
   }

   return 0;
}
