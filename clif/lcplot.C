#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <libsite.H>
#include "libclif.H"
#include "supermongo.h"

extern "C" {
   void sm_defvar( char *, char * );
   void sm_ltype( float );
   void sm_expand( double );
   void sm_ticksize( float, float, float, float );
   void sm_lweight( int );
   void sm_angle( double );
}

int LcDataLimits( CommonLc &lc, float &xlimlo, float &xlimhi,
                  float &ylimlo, float &ylimhi, float &elimlo, float &elimhi )
{
   for( int n = 0; n < lc.nobs; ++n ) {
      xlimlo = lc.date[n] < xlimlo ? lc.date[n] : xlimlo; 
      xlimhi = lc.date[n] > xlimhi ? lc.date[n] : xlimhi;
      ylimlo = lc.value[n] < ylimlo ? lc.value[n] : ylimlo; 
      ylimhi = lc.value[n] > ylimhi ? lc.value[n] : ylimhi;
      elimlo = lc.error[n] < elimlo ? lc.error[n] : elimlo; 
      elimhi = lc.error[n] > elimhi ? lc.error[n] : elimhi;
   }
   return SiteErr::errNone;
}
 
int LcDataLimits( CommonLc &lc, float &xlimlo, float &xlimhi,
                  float &ylimlo, float &ylimhi )
{
   static float elimlo, elimhi;
   return LcDataLimits( lc, xlimlo, xlimhi, ylimlo, ylimhi, elimlo, elimhi );
}
 
int LcDataLimits( LcSet &s, float &xlimlo, float &xlimhi,
                  float &ylimlo, float &ylimhi, float &elimlo, float &elimhi )
{
   xlimlo = ylimlo = elimlo = Infinity;
   xlimhi = ylimhi = elimhi = -Infinity;
   for( int p = 0; p < s.count( ); ++p ) {
      LcDataLimits( *s[p], xlimlo, xlimhi, ylimlo, ylimhi, elimlo, elimhi );
   }
   return SiteErr::errNone;
}
 
int LcDataLimits( LcSet &s, float &xlimlo, float &xlimhi,
                  float &ylimlo, float &ylimhi )
{
   static float elimlo, elimhi;
   return LcDataLimits( s, xlimlo, xlimhi, ylimlo, ylimhi, elimlo, elimhi );
}

LcPlot::LcPlot( void ) {
   ok = maxn = nobs = maxp = 0;
   xb = yb = eb = 0;
   *devstr = *xlabel = *ylabel = *plabel = 0;
   xlimlo = ylimlo = Infinity;
   xlimhi = ylimhi = -Infinity;
   ylim[0] = ylim[1] = 0;
}

LcPlot::~LcPlot( void ) {
   delete [] xb;
   delete [] yb;
   delete [] eb;
   delete [] ylim[0];
   delete [] ylim[1];
}

int LcPlot::resize( int size ) {
   if( size > maxn ) {
      delete [] xb;
      delete [] yb;
      delete [] eb;
      if( !(xb = new float[size]) ||
	  !(yb = new float[size]) ||
          !(eb = new float[size]) )
	 return SiteErr::errMemory;
   }
   return SiteErr::errNone;
}

int LcPlot::loadlc( CommonLc &lc ) {
   if( lc.nobs == 0 )
      return SiteErr::errMissing;
   if( resize( lc.nobs ) )
      return SiteErr::errMemory;
   nobs = lc.nobs;
   int n;
   for( n = 0; n < nobs; ++n ) {
      xb[n] = lc.date[n];
      yb[n] = lc.value[n];
      eb[n] = lc.error[n];
   }
   return SiteErr::errNone;
}

int LcPlot::config( LcSet &s ) {

   *devstr = *xlabel = *ylabel = *plabel = 0;
   sprintf( devstr, "x11 -bg white -fg black -geom 600x600+0+0" );

   if( s.count( ) == 0 )
      return SiteErr::errMissing; 

   /*--- plot limits ---*/
   xlimlo = ylimlo = Infinity;		// redundant with LcDataLimits
   xlimhi = ylimhi = -Infinity;		// redundant with LcDataLimits
   LcDataLimits( s, xlimlo, xlimhi, ylimlo, ylimhi );
   float dy = (ylimhi - ylimlo);
   ylimlo -= 0.2 * dy;			// expand limits
   ylimhi += 0.2 * dy;			// expand limits

   /*--- set y limits for individual lightcurves ---*/
   /*--- only data lightcurve are counted and sized ---*/
   if( s.count( ) > maxp ) {
      delete [] ylim[0];
      delete [] ylim[1];
      if( !(ylim[0] = new float[s.count( )]) ||
	  !(ylim[1] = new float[s.count( )]) ) 
	 return SiteErr::errMemory;
   }
   int d, n;
   float xd0, xd1;
   for( d = n = 0; n < s.count( ); ++n ) {
      if( (*s[n]).phtpkg == CommonLc::Theory ) 
	 continue;
      ylim[1][d] = - (ylim[0][d] = Infinity);
      LcDataLimits( *s[n], xd0, xd1, ylim[0][d], ylim[1][d]);
      dy =  ylim[1][d] - ylim[0][d];
      ylim[0][d] -= dy * 0.2;
      ylim[1][d] += dy * 0.2;
      ++d;			// data curve number
   }
   
   /*--- plot labels ---*/
   sprintf( plabel, "%s", s[0]->title );
   sprintf( xlabel, "JD - %.6f", s[0]->refJD );
   if( s[0]->units == CommonLc::Magnitude ) 
      sprintf( ylabel, "Magnitudes" );
   else
      sprintf( ylabel, "Linear Units" );
   ynot = 0.9;				// notation expand paramter
   xnot = 0.9;				// notation expand paramter

   /*--- point type ---*/
   ptype[0] = 8; ptype[1] = 3;		// ptype 8 3   (solid octagons)

   ok = 1;
   return SiteErr::errNone; 
}

int LcPlot::draw( LcSet &s ) {

   int nt = 0, nd = 0;				// number of theory & data Lc's
   int wn = 0, t, n, d;
   int setlim = 0;
   char msg[256];
   CommonLc *lcp, *tlc[64];
   float xlab, ylab;
   float dx, dy;
   float pxlimlo, pxlimhi, pylimlo, pylimhi;

   sm_defvar("TeX_strings", "1");
   if( !ok )
      return SiteErr::errProgram;

   sm_device( devstr );
   sm_graphics( );
   sm_erase( );
   sm_ptype( ptype, 2 );
   sm_device( devstr );
   //sm_expand( 1.05 );

   if( s.count( ) == 0 )
      return 0;

   /*--- search for theory curves and hang onto pointers ---*/
   for( n = 0; n < s.count( ); ++n )
      if( (*s[n]).phtpkg == CommonLc::Theory ) {
	 if( (*s[n]).nobs > 0 ) 
	    tlc[nt++] = s[n];
      } else
	 ++nd;


   //fprintf(stderr,"LcPlot::draw %d data curves, %d theory curves\n",nd,nt);
   //fprintf(stderr, "LcPlot::draw devstr: %s\n", devstr );
   //fprintf(stderr, "LcPlot::draw plabel: %s\n", plabel );
   //fprintf(stderr, "LcPlot::draw xlabel: %s\n", xlabel );
   //fprintf(stderr, "LcPlot::draw ylabel: %s\n", ylabel );
   //fprintf(stderr, "LcPlot::draw xlimits %.4g -- %.4g\n", xlimlo, xlimhi );
   //fprintf(stderr, "LcPlot::draw ylimits %.4g -- %.4g\n", ylimlo, ylimhi );
   dx = xlimhi - xlimlo;	// t limits fixed from here on
   dy = ylimhi - ylimlo;	// not used anymore
   wn = nd;
   for( d = t = n = 0; n < s.count( ); ++n ) {
      
      if( (*s[n]).phtpkg == CommonLc::Theory )
	 continue;
      lcp = s[n];
      sm_window( 1, -nd, 1, wn, 0, 0 );
 
      if( lcp->units == CommonLc::Magnitude )
	 sm_limits( xlimlo, xlimhi, ylim[1][d], ylim[0][d]);
      else
	 sm_limits( xlimlo, xlimhi, ylim[0][d], ylim[1][d]);
      dy = ylim[1][d] - ylim[0][d];
      
      if( nt && t != nt ) {
	 loadlc( *tlc[t] );
	 sm_conn( xb, yb, nobs );
	 ++t;
      }

      /*--- color data ---*/
      float w = Passband::info( (*s[n]).filter )->midpoint;
      if( w == 5200 )			// macho blue is blue regardless
	 sm_ctype( "blue" );
      else if( w > 8000 )		// I band
	 sm_ctype( "magenta" );
      else if( w > 6000 )		// R band
	 sm_ctype( "red" );
      else if( w > 5000 )		// V band
	 sm_ctype( "green" );
      else if( w > 3500 )		// B band
	 sm_ctype( "blue" );
      else
	 sm_ctype( "default" );

      loadlc( *lcp );
      sm_errorbar( xb, yb, eb, 2, nobs );
      sm_errorbar( xb, yb, eb, 4, nobs );
      sm_points( xb, yb, nobs );
      sm_ctype( "default" );
      
      xlab = xlimlo + dx * xnot;
      if( lcp->units == CommonLc::Magnitude )
	 ylab = ylim[0][d] + dy * (1-ynot);
      else
	 ylab = ylim[0][d] + dy * ynot;
      sm_relocate( xlab, ylab );
      sm_putlabel( 4, (*Passband::info( lcp->filter ) ).name );
      sm_box( 0, 2, 0, 0 );
      ++d;
      --wn;
   }
   sm_window( 1, 1, 1, 1, 0, 0 );
   sm_limits( xlimlo, xlimhi, 0, 1 );
   sm_box( 1, 3, 0, 3 );
   sm_limits( 0, 1, 0, 1 );
   xlab = 0.5;
   ylab = 1.03;
   sm_relocate( xlab, ylab );
   if( *plabel )
      sm_putlabel( 5, plabel );
   if( *xlabel )
      sm_xlabel( xlabel );
   xlab = -0.11;
   ylab = 0.5;
   sm_angle( 90.0 );
   sm_relocate( xlab, ylab );
   if( *ylabel )
      sm_putlabel( 5, ylabel );
   sm_angle( 0.0 );
   sm_gflush( );
   sm_hardcopy( );
   return SiteErr::errNone;
}

int LcPlot::print( LcSet &s, const char *dev, const char *filename ) {
   if( filename ) sprintf( devstr, "%s %s", dev, filename );
   else sprintf( devstr, "%s", dev );
   draw( s );

   /*--- reset to screen default ---*/
   sprintf( devstr, "x11 -bg white -fg black -geom 600x600+0+0" );
   return SiteErr::errNone;
}

