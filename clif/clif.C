/* written by John Doug Reynolds, April 1998 */

#include "libclif.H"

#include <libsite.H>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>

ClifParameters clifpar;			// global set of defaults

/* --- CommonLc methods ---------------------------------------------------- */

CommonLc::CommonLc( void )
{
   maxobs = 0;
   title = 0;
   date = value = error = 0;
   obsid = 0;
   erase();
}

CommonLc::CommonLc( const char *TITLE, int NOBS )
{
   maxobs = 0;
   title = 0;
   date = value = error = 0;
   obsid = 0;
   erase();
   title = duplicate( TITLE );
   resize( NOBS );
}

CommonLc::~CommonLc( void )
{
   delete [] title;
   delete [] date;
   delete [] value;
   delete [] error;
   delete [] obsid;
}

void	CommonLc::copy( const CommonLc &lc, int header_only )
{
   if ( title ) delete [] title;

   title	= duplicate( lc.title );
   ra		= lc.ra;
   dec		= lc.dec;
   refJD	= lc.refJD;
   mag0		= lc.mag0;
   norm1	= lc.norm1;
   norm0	= lc.norm0;
   units	= lc.units;
   filter	= lc.filter;
   phtpkg	= lc.phtpkg;
   id		= lc.id;

   if ( !header_only && !resize(lc.nobs) && nobs > 0 ) {
      memcpy( date, lc.date, nobs * sizeof(double) );
      memcpy( value, lc.value, nobs * sizeof(double) );
      memcpy( error, lc.error, nobs * sizeof(double) );
      memcpy( obsid, lc.obsid, nobs * sizeof(int) );
   }
   else nobs = 0;
}

void	CommonLc::erase( void )
{
   if ( title ) delete [] title;
   title = 0;
   ra = dec = NaN;
   refJD = mag0 = norm1 = norm0 = 0;
   units = 0;
   filter = Passband::Invalid;
   phtpkg = Undefined;
   id = -1;
   nobs = 0;
}

int	CommonLc::resize( int NOBS )
{
   if ( NOBS < 0 ) return 1;

   if ( NOBS > maxobs ) {
      double *ftemp;
      if ( !(ftemp = new double [NOBS]) ) return 1;
      if ( date ) {
	 memcpy( (char*) ftemp, (char*) date, (int) nobs * sizeof(double) );
	 delete [] date;
      }
      date = ftemp;

      if ( !(ftemp = new double [NOBS]) ) return 1;
      if ( value ) {
	 memcpy( (char*) ftemp, (char*) value, (int) nobs * sizeof(double) );
	 delete [] value;
      }
      value = ftemp;

      if ( !(ftemp = new double [NOBS]) ) return 1;
      if ( error ) {
	 memcpy( (char*) ftemp, (char*) error, (int) nobs * sizeof(double) );
	 delete [] error;
      }
      error = ftemp;

      int *ltemp;
      if ( !(ltemp = new int [NOBS]) ) return 1;
      if ( obsid ) {
	 memcpy( (char*) ltemp, (char*) obsid, (int) nobs * sizeof(int) );
	 delete [] obsid;
      }
      obsid = ltemp;

      maxobs = NOBS;
   }

   nobs = NOBS;
   return 0;
}

int	CommonLc::convert( Units u )
{
   if ( units == u ) return 0;

   // the error conversions here are approximate and need more thought

   switch ( u ) {
   default: return 1;
   case Magnitude:
      {
	 const double C = 2.5 / log( 10.0 );
	 //const double N = norm0 / norm1;
	 for ( int o = 0; o < nobs; ++o ) {
	    error[o] = C * error[o] / value[o];
	    value[o] = -2.5 * log10( value[o] ) + norm1;
	 }
	 //norm0 = -2.5 * log10( norm1 );
	 //norm1 = 1;
	 break;
      }
   case Linear:
      {
	 const double C = log( 10.0 ) / 2.5;
	 //const double N = -0.4 * norm1;
	 for ( int o = 0; o < nobs; ++o ) {
	    value[o] = pow( 10.0, -0.4 * (value[o] - norm1) );
	    error[o] = C * error[o] * value[o];
	 }
	 //norm1 = pow( 10.0, -0.4 * norm0 );
	 //norm0 = 0;
	 break;
      }
   }

   units = u;
   return 0;
}   

int	CommonLc::renorm( double n1, double n0 )
{
   if ( n1 == 0 ) return 1;
   if ( units == Magnitude  &&  n1 != 1 ) return 1;

   int o;
   const double invn1 = 1 / n1;
   for ( o = 0; o < nobs; ++o ) value[o] = invn1 * ( value[o] - n0 );
   if ( units == Linear )
      for ( o = 0; o < nobs; ++o ) error[o] = invn1 * error[o];

   norm0 += norm1 * n0;
   norm1 *= n1;
   return 0;
}

int	CommonLc::print( FILE *fi, int start, int count ) const
{
   if ( start < 0 ) start = 0;
   int max = count < 0 ? nobs : start + count;
   if ( max > nobs ) max = nobs;

   char *ustr, *pkg;
   switch ( phtpkg ) {
   default:		pkg = "unknown"; break;
   case Undefined:	pkg = "undefined"; break;
   case Theory:		pkg = "theory"; break;
   case pSodophot:	pkg = "production SoDOPHOT"; break;
   case oSodophot:	pkg = "offline SoDOPHOT"; break;
   case Sodophot:	pkg = "SoDOPHOT"; break;
   case Daophot:	pkg = "DaoPhot"; break;
   case Allframe:	pkg = "Allframe"; break;
   case DIPhot:		pkg = "Image Difference"; break;
   }
   switch ( units ) {
   default:		ustr = "unknown"; break;
   case Magnitude:	ustr = "magnitude"; break;
   case Linear:		ustr = "linear"; break;
   }
   const Passband::Info *info = Passband::info( filter );
   const double mas_rad = 648e6 / M_PI;

   fprintf( fi, "! %s\n#\n", title );
   fprintf( fi, "# %s %s (J2000.0)\n"
	    , macho_ra_string( isnan(ra) ? 0 : int(0.5+ra*mas_rad), ':' )
	    , macho_dec_string( isnan(dec) ? 0 : int(0.5+dec*mas_rad), ':' ) );
   fprintf( fi, "# photometry in %s reduced by %s\n"
	    , info ? info->name : "no band", pkg );
   fprintf( fi, "# data in %s units\n", ustr );
   fprintf( fi, "# ADU = %.3f * value + %.3f\n", norm1, norm0 );
   fprintf( fi, "# magnitude zero-point: %.3f\n", mag0 );
   fprintf( fi, "# reference Julian date: %.15g\n", refJD );
   fprintf( fi, "# %d observations\n", count ? max-start : nobs );
   fprintf( fi, "#\n!%14s%15s%15s%15s\n", "date", "value", "error", "obsid" );

   for ( int o = start; o < max; ++o )
      if ( 0 > fprintf( fi, "%15.8f%#15.8g%#15.8g%15ld\n"
			, date[o], value[o], error[o], obsid[o] ) ) return 1;

   return 0;
}

/* --- Passband methods ---------------------------------------------------- */

static Passband::Info passbandinfo[] = {
   { 0,    0,	 "Invalid",	 "INVALID",	Passband::Invalid },
   { 7000, 2200, "Standard R",	 "STND.R",	Passband::R_Std },
   { 5500, 890,  "Standard V",	 "STND.V",	Passband::V_Std },
   { 4400, 980,  "Standard B",	 "STND.B",	Passband::B_Std },
   { 6850, 1900, "MACHO Red",	 "MACHO.R",	Passband::R_MACHO },
   { 5200, 1400, "MACHO Blue",	 "MACHO.V",	Passband::V_MACHO },
   { 9000, 2400, "CTIO 0.9m I",	 "CTIO.I",	Passband::I_CTIO },
   { 7000, 2200, "CTIO 0.9m R",	 "CTIO.R",	Passband::R_CTIO },
   { 5500, 890,  "CTIO 0.9m V",	 "CTIO.V",	Passband::V_CTIO },
   { 4400, 980,  "CTIO 0.9m B",	 "CTIO.B",	Passband::B_CTIO },
   { 3650, 680,  "CTIO 0.9m U",	 "CTIO.U",	Passband::U_CTIO },
   { 9000, 2400, "WISE 1m I",	 "WISE.I",	Passband::I_WISE },
   { 7000, 2200, "WISE 1m R",	 "WISE.R",	Passband::R_WISE },
   { 5500, 890,  "WISE 1m V",	 "WISE.V",	Passband::V_WISE },
   { 4400, 980,  "WISE 1m B",	 "WISE.B",	Passband::B_WISE },
   { 3650, 680,  "WISE 1m U",	 "WISE.U",	Passband::U_WISE },
   { 9000, 2400, "MJUO 0.61m I", "MJUO.I",	Passband::I_MJUO },
   { 7000, 2200, "MJUO 0.61m R", "MJUO.R",	Passband::R_MJUO },
   { 5500, 890,  "MJUO 0.61m V", "MJUO.V",	Passband::V_MJUO },
   { 4400, 980,  "MJUO 0.61m B", "MJUO.B",	Passband::B_MJUO },
   { 3650, 680,  "MJUO 0.61m U", "MJUO.U",	Passband::U_MJUO },
   { 9000, 2400, "UTSO 0.61m I", "UTSO.I",	Passband::I_UTSO },
   { 7000, 2200, "UTSO 0.61m R", "UTSO.R",	Passband::R_UTSO },
   { 5500, 890,  "UTSO 0.61m V", "UTSO.V",	Passband::V_UTSO },
   { 4400, 980,  "UTSO 0.61m B", "UTSO.B",	Passband::B_UTSO },
   { 3650, 680,  "UTSO 0.61m U", "UTSO.U",	Passband::U_UTSO },
   { 9000, 2400, "MSO 0.8m I",	 "MSO30.I",	Passband::I_MSO30 },
   { 7000, 2200, "MSO 0.8m R",	 "MSO30.R",	Passband::R_MSO30 },
   { 5500, 890,  "MSO 0.8m V",	 "MSO30.V",	Passband::V_MSO30 },
   { 4400, 980,  "MSO 0.8m B",	 "MSO30.B",	Passband::B_MSO30 },
   { 3650, 680,  "MSO 0.8m U",	 "MSO30.U",	Passband::U_MSO30 },
   { 9000, 2400, "MSO 1.9m I",	 "MSO74.I",	Passband::I_MSO74 },
   { 7000, 2200, "MSO 1.9m R",	 "MSO74.R",	Passband::R_MSO74 },
   { 5500, 890,  "MSO 1.9m V",	 "MSO74.V",	Passband::V_MSO74 },
   { 4400, 980,  "MSO 1.9m B",	 "MSO74.B",	Passband::B_MSO74 },
   { 3650, 680,  "MSO 1.9m U",	 "MSO74.U",	Passband::U_MSO74 },
   { 9000, 2400, "APO 3.5m I",	 "APO.I",	Passband::I_APO },
   { 7000, 2200, "APO 3.5m R",	 "APO.R",	Passband::R_APO },
   { 5500, 890,  "APO 3.5m V",	 "APO.V",	Passband::V_APO },
   { 4400, 980,  "APO 3.5m B",	 "APO.B",	Passband::B_APO },
   { 3650, 680,  "APO 3.5m U",	 "APO.U",	Passband::U_APO },
   { 9000, 2400, "CTIO 1.5m I",	 "CTIO15.I",	Passband::I_CTIO15 },
   { 7000, 2200, "CTIO 1.5m R",	 "CTIO15.R",	Passband::R_CTIO15 },
   { 5500, 890,  "CTIO 1.5m V",	 "CTIO15.V",	Passband::V_CTIO15 },
   { 4400, 980,  "CTIO 1.5m B",	 "CTIO15.B",	Passband::B_CTIO15 },
   { 3650, 680,  "CTIO 1.5m U",	 "CTIO15.U",	Passband::U_CTIO15 },
   { 6467, 1509, "VATT 1.8m RJT","VATT.RJT",	Passband::RJT_VATT }
};

const Passband::Info*	Passband::info( const char *tag )
{
   if ( strlen(tag) < 16 ) {
      int i;
      char TAG[16];
      const int max = sizeof( passbandinfo ) / sizeof( Passband::Info );
      for ( i = 0; i<16 && (TAG[i]=toupper(tag[i])); ++i );
      for ( i = 0; i < max; ++i )
	 if ( strcmp(passbandinfo[i].tag,TAG) == 0 ) return &passbandinfo[i];
   }
   return &passbandinfo[0];
}

const Passband::Info*	Passband::info( short id )
{
   const int max = sizeof( passbandinfo ) / sizeof( Passband::Info );
   for ( int i = 0; i < max; ++i )
      if ( passbandinfo[i].id == id ) return &passbandinfo[i];
   return &passbandinfo[0];
}

/* --- LcSet methods ------------------------------------------------------- */

LcSet::LcSet( void )
{
   numrip = maxrip = 0;
   rip = 0;
   num = max = 0;
   lc = 0;
}

LcSet::~LcSet( void )
{
   int i;
   if ( lc ) for ( i = 0; i < num; ++i ) delete lc[i];
   delete [] lc;
   if ( rip ) for ( i = 0; i < numrip; ++i ) delete rip[i];
   delete [] rip;
}

CommonLc*	LcSet::operator[] ( int index ) const
{
   return index >= 0 ? index < num ? lc[index] : 0 : 0;
}

int	LcSet::append( CommonLc *clc )
{
   for ( int i = 0; i < num; ++i ) if ( lc[i] == clc ) return i;

   if ( num == max ) {
      CommonLc **temp = new CommonLc* [max+10];
      if ( !temp ) return -1;
      max += 10;
      if ( lc ) {
	 memcpy( (char*) temp, (char*) lc, num * sizeof( CommonLc* ) );
	 delete [] lc;
      }
      lc = temp;
   }

   lc[num++] = clc;
   return num;
}   

CommonLc*	LcSet::newlc( void )
{
   CommonLc *clc = numrip ? rip[--numrip] : new CommonLc;
   if ( !clc ) return 0;
   clc->erase();
   return clc;
}

int	LcSet::unhook( const CommonLc *clc )
{
   int i;
   for ( i = 0; i < num; ++i ) if ( lc[i] == clc ) break;
   if ( i < num ) {
      for (; i+1 < num; ++i ) lc[i] = lc[i+1];
      lc[--num] = 0;
      return 0;
   }
   return 1;
}

int	LcSet::recycle( CommonLc *clc )
{
   if ( !clc  ||  unhook(clc) ) return 1;

   if ( numrip == maxrip ) {
      CommonLc **temp = new CommonLc* [maxrip+10];
      if ( !temp ) return 1;
      maxrip += 10;
      if ( rip ) {
	 memcpy( (char*) temp, (char*) rip, numrip * sizeof( CommonLc* ) );
	 delete [] rip;
      }
      rip = temp;
   }

   rip[numrip++] = clc;
   return 0;
}

int	LcSet::recycle( void )
{
   while ( num )
      if ( recycle( lc[0] ) ) return 1;
}

/* --- ClifParameters methods ---------------------------------------------- */

ClifParameters::ClifParameters( void )
{
   num = max = 0;
   item = 0;
   filename = 0;
}

ClifParameters::~ClifParameters( void )
{
   reset();
   delete [] item;
}

void	ClifParameters::reset( void )
{
   for ( int i = 0; i < num; ++i ) {
      Pair &p = item[i];
      delete [] p.name;
      delete [] p.value;
      p.name = p.value = NULL;
   }
   num = 0;
   delete [] filename;
   filename = 0;
}

// binary search
int	ClifParameters::at( const char *name ) const
{
   if ( !name ) return -1;

   int imin = 0;
   int imax = num - 1;

   while ( imax >= imin ) {
      int i = (imax + imin) / 2;
      int diff = strcmp( item[i].name, name );
      if ( diff == 0 ) return i;
      else if ( diff < 0 ) imin = i + 1;
      else imax = i - 1;
   }

   return imin;
}

int	ClifParameters::add( const char *n, const char *v )
{
   char *value = new char [ strlen(v) + 1 ];
   if ( !value ) {
      fprintf( stderr, "ClifParameters: out of memory\n" );
      return 1;
   }
   strcpy( value, v );

   int index = at(n);
   if ( index < num  &&  strcmp(n,item[index].name) == 0 ) {
      delete [] item[index].value;
      item[index].value = value;
      //fprintf( stderr, "ClifParameters reassigning %s %s\n", n, v );  
      return 0;
   }

   if ( num == max ) {
      Pair *temp = new Pair [max+256];
      if ( !temp ) {
	 fprintf( stderr, "ClifParameters: out of memory\n" );
	 return 1;
      }
      max += 256;
      if ( item ) {
	 memcpy( (char*) temp, (char*) item, sizeof(Pair) );
	 delete [] item;
      }
      item = temp;
   }

   char *name = new char [ strlen(n) + 1 ];
   if ( !name ) {
      //fprintf( stderr, "ClifParameters: out of memory\n" );
      return 1;
   }
   strcpy( name, n );

   for ( int i = num - 1; i >= index; --i ) item[i+1] = item[i];
   item[index].name = name;
   item[index].value = value;
   ++num;
   //fprintf( stderr, "ClifParameters assigning %s %s\n", n, v );  

   return 0;
}

const ClifParameters::Pair*   ClifParameters::operator[] ( int index ) const
{
   return index < num ? index >= 0 ? &item[index] : 0 : 0;
}

const char*	ClifParameters::operator[] ( const char *name ) const
{
   int i = at( name );
   if ( i < 0  ||  i >= num  ||  strcmp(item[i].name,name) != 0 ) return NULL;
   return item[i].value;
}

int	ClifParameters::get( const char *name, int def ) const
{
   int param;
   const char *value = (*this)[name];
   if ( value  &&  sscanf( value, "%i", &param ) == 1 ) return param;
   return def;
}

double	ClifParameters::get( const char *name, double def ) const
{
   double param;
   const char *value = (*this)[name];
   if ( value  &&  sscanf( value, "%lf", &param ) == 1 ) return param;
   return def;
}

int	ClifParameters::load( const char *fname )
{
   if ( !fname  &&  !(fname=getenv("CLIFPAR")) ) {
      char *name = ".clif-parameters";
      char *HOME = getenv( "HOME" );
      if ( !HOME ) {
	 fprintf(stderr,"ClifParameters: missing environment variable HOME\n");
	 return 1;
      }
      filename = new char [ strlen(HOME) + strlen(name) + 1 ];
      if ( !filename ) {
	 fprintf( stderr, "ClifParameters: out of memory\n" );
	 return 1;
      }
      sprintf( filename, "%s/%s", HOME, name );
      struct stat sbuf;
      if ( stat(filename,&sbuf)  &&  errno == ENOENT ) return 0;
   }
   else {
      if ( !*fname ) return 0;
      filename = new char [ strlen(fname) + 1 ];
      if ( !filename ) {
	 fprintf( stderr, "ClifParameters: out of memory\n" );
	 return 1;
      }
      strcpy( filename, fname );
   }

   FILE *fi = fopen( filename, "r" );
   if ( !fi ) {
      fprintf( stderr, "ClifParameters: could not open: %s\n", filename );
      return 1;
   } 
   //fprintf( stderr, "ClifParameters reading %s\n", filename );

   int err = 0;
   char buf[256], name[256];
   for ( int line = 1; !err && fgets( buf, 256, fi ); ++line ) {
      // strip off trailing junk
      char *ptr = strchr( buf, '#' );
      if ( ptr ) *ptr = '\0';
      ptr = buf + strlen(buf) - 1;
      while ( ptr >= buf  &&  isspace(*ptr) ) *ptr-- = '\0';

      // read parameter name and strip off leading junk
      int n;
      if ( sscanf( buf, "%s%n", name, &n ) < 1 ) continue;
      for ( ptr = buf+n; isspace(*ptr); ++ptr );

      // check for a value before adding the pair
      if ( !*ptr ) {
	 fprintf(stderr,"ClifParameters: %s:%d missing value\n",filename,line);
	 continue;
      }
      err = add( name, ptr );
   }

   fclose( fi );
   return err;
}
