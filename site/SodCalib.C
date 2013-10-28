/* written by John Doug Reynolds, February 1998 */

#include "Macho.h"
#include "libsite.H"
#include <stdlib.h>
#include <string.h>
#include <errno.h>

SodCalibrator::SodCalibrator( void )
{
   rwcid = 255;
   field = obsid = npsf = 0;
   airmass = a0 = a1 = b0 = b1 = co = mc = NaN;
   bje = bjw = bjo = V0 = V1 = R0 = R1 = NaN;
}

SodCalibDb::SodCalibDb( void )
{
   field = 0;
   memset( chip, 0, sizeof(chip) );
   memset( chunk, 0, sizeof(chunk) );
}

int	SodCalibDb::read( int field_, const char *lcdir, const char *active )
{
   if ( field_ <= 0 || (!lcdir && !active) ) return SiteErr::errParams;

   int err;
   char fname[1024];

   if ( lcdir ) {
      sprintf(fname,"%s/Calib/SodCalib_%03d",lcdir,field_);
      if ( ! (err = read(fname)) ) return 0;
      if ( err != -ENOENT ) return err;
   }

   if ( active ) {
      sprintf(fname,"%s/StarByField/F_%d/SodCalib_%03d",active,field_,field_);
      if ( ! (err = read(fname)) ) return 0;
      if ( err != -ENOENT ) return err;
   }

   return SiteErr::errMissing;
}

int	SodCalibDb::read( const char *filename )
{
   int i, err = 0;
   char scratch[1024];
   FILE *fi = fopen(filename,"r");
   if ( ! fi ) return -errno;

   if ( !fgets(scratch,1024,fi)			||
	!fgets(scratch,1024,fi)			||
	sscanf(scratch,"field %d",&field) != 1	||
	!fgets(scratch,1024,fi)
	)
      err = SiteErr::errFormat;

   for ( i = 0; !err && i < 4; ++i )
      if ( fgets(scratch,1024,fi) ) {
	 int index = atoi(scratch);
	 if ( index >= 0 && index < 4 ) {
	    Chip &s = chip[index];
	    if ( 6 != sscanf(scratch,"%*d%*d%d%f%f%f%f%f"
			     ,&s.obsid,&s.airmass,&s.a0,&s.a1,&s.b0,&s.b1) )
	       err = SiteErr::errFormat;
	 }
	 else err = SiteErr::errFormat;
      }
      else err = SiteErr::errFormat;

   if ( !err && !fgets(scratch,1024,fi) ) err = SiteErr::errFormat;

   for ( i = 0; !err && i < 64; ++i )
      if ( fgets(scratch,1024,fi) ) {
	 int rwcid = atoi(scratch);
	 if ( (rwcid >= 0 && rwcid < 56) || (rwcid >= 120 && rwcid < 128) ) {
	    Chunk &s = chunk[(rwcid+8)&63];
	    s.rwcid = rwcid;
	    if ( 4 != sscanf(scratch,"%*d%f%f%f%d",&s.co,&s.bj,&s.mc,&s.npsf) )
	       err = SiteErr::errFormat;
	 }
	 else err = SiteErr::errFormat;
      }
      else err = SiteErr::errFormat;

   fclose(fi);
   if ( err ) field = 0;
   return err;
}

int	SodCalibDb::get( class SodCalibrator &calib, int rwcid ) const
{
   if ( ! field ) return SiteErr::errCorrupt;
   if ( rwcid < 0 || (rwcid > 55 && rwcid < 120) || rwcid > 127 )
      return SiteErr::errParams;

   int index = (rwcid + 8) & 63;
   int chipid = index >> 4;

   const Chunk &ck = chunk[index];
   const Chip &cp = chip[chipid];

   calib.field = field;
   calib.rwcid = rwcid;
   calib.obsid = cp.obsid;
   calib.npsf = ck.npsf;
   calib.airmass = cp.airmass;
   calib.a0 = cp.a0;
   calib.a1 = cp.a1;
   calib.b0 = cp.b0;
   calib.b1 = cp.b1;
   calib.co = ck.co;
   calib.mc = ck.mc;

   int swap = (field < 0) ? chipid == 2 : chipid == 0;
   calib.bje = swap ? -ck.bj : 0;
   calib.bjw = swap ? 0 : ck.bj;

   calib.bjo = (ck.mc < 0.35) ? -0.05 : ck.mc - 0.4;

   calib.V0 = cp.a0 + ck.co;
   calib.V1 = cp.a1 + 0.022 * cp.airmass;
   calib.R0 = cp.b0 + ck.co;
   calib.R1 = cp.b1 + 0.004 * cp.airmass;

   return 0;
}

int	SodCalibDb::print( FILE *fi ) const
{
   time_t clock = time(0);
   struct tm *tmp = localtime( &clock );
   fprintf( fi, "# SoDophot calibration database version: %02d%02d%02d\n"
	    , tmp->tm_year, tmp->tm_mon+1, tmp->tm_mday );
   fprintf( fi, "field %-5d# %s-style template\n"
	    , field, field < 0 ? "east" : "west" );

   int i, zpc[4] = {3,19,31,51};
   char *rwcids[4] = {"0-7,120-127","8-23","24-39","40-55"};
   fprintf( fi, "#%6s%6s%7s%9s%8s%8s%8s%8s%5s  %s\n"
	    , "rchip", "bchip", "obsid", "airmass"
	    , "a0", "a1", "b0", "b1", "zpc", "rwcids" );
   for ( i = 0; i < 4; ++i ) {
      const SodCalibDb::Chip &s = chip[i];
      fprintf( fi, "%7d%6d%7d%9.5f%8.3f%8.4f%8.3f%8.4f%5d  %s\n"
	       , i, 7-i, s.obsid, s.airmass
	       , s.a0, s.a1, s.b0, s.b1, zpc[i], rwcids[i] );
   }

   fprintf( fi, "#%6s%10s%10s%10s%8s\n"
	    , "rwcid", "co", "bj", "mc_psf", "no_psf" );
   for ( i = 0; i < 64; ++i ) {
      const SodCalibDb::Chunk &s = chunk[i];
      fprintf( fi, "%7d%10.3f%10.3f%10.3f%8d\n"
	       , s.rwcid, s.co, s.bj, s.mc, s.npsf );
   }

   return 0;
}
