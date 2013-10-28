#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <libsite.H>
#include <libclif.H>
#include "lcfit.H"


int main( int argc, char *argv[] )
{
   if( argc < 2 ) {
      printf("\nUsage:\n\t%s photfile\n", argv[0] );
      exit(1);
   }

   FILE *fpin;
   CommonLc *lcp;
   LcSet lcs, hlcs;
   PeakFinder pf( clifpar );
   LcMoment mom;
   char ibuf[1024];		// input char buffer

   int n, r;
   if( (fpin = fopen( argv[1], "r" ) ) == 0 ) {
      printf("Could not read %s, exiting\n", argv[n] );
      return 1;
   }
   while( fgets( ibuf, 1024, fpin ) ) ++r;	// count lines
   rewind( fpin );
   lcs.append( lcp = lcs.newlc( ) );		// append lc to set
   lcp->resize( r );
   while( fgets( ibuf, 1024, fpin ) ) {
      if( *ibuf == '#' ) continue;		// skip comments
      sscanf( ibuf, "%lf %lf %lf",
	      lcp->date+r, lcp->value+r, lcp->error+r );
      ++r;
   }
   lcp->nobs = r;
   lcp->filter = Passband::RJT_VATT;
   lcp->phtpkg = CommonLc::DIPhot;
   lcp->units = CommonLc::Magnitude;
   lcp->title = strdup( argv[1] );
   lcp->convert( CommonLc::Linear );	// and convert to linear units

   if( !lcs.count( ) ) {
      printf("No data: exiting\n" );
      return 1;
   }

   if( lcs[0]->nobs < 5 ) {
      printf("Not enough data (%d points): exiting\n", lcs[0]->nobs );
      return 1;
   }
   printf("Loaded %d points from %s\n", lcs[0]->nobs,lcs[0]->title );

   mom.eval( *lcs[0] );
   lcs[0]->renorm( mom.average( ), 0.0 );	// renormalize to average

   pf.eval( lcs );	/*--- run PeakFinder on LcSet ---*/
   MicroBlendFit fmicro( clifpar );
   fmicro.init( lcs, pf );

   /*--- Append theory Lc's for later use ---*/
   CommonLc *tlc;
   lcs.append( tlc = new CommonLc );
   lc_hdr_cpy( *lcs[0], *tlc );
   tlc->units = CommonLc::Linear;
   tlc->phtpkg = CommonLc::Theory;

   /*--- Show raw data for parameter guessing ---*/
   LcPlot plot;
   plot.config( lcs );
   plot.draw( lcs );
   printf("\t> " );   gets( ibuf );	// pause for keystroke

   fmicro.iopar( 1 );				// query for initial values
   fmicro.ifit( 0, 0, 1 );			// do fit
   fmicro.iopar( );				// print results

   /*--- make theory curve ---*/
   lcs[1]->resize( 1001 );
   double dt = (plot.xlimhi - plot.xlimlo) / 1000.0;
   for( n = 0; n < 1001; ++n )
      lcs[1]->date[n] = plot.xlimlo + dt*n;		
   fmicro.FitLc( *lcs[1], 0 );

   plot.config( lcs );
   plot.draw( lcs );
   printf("Exit ? " );
   gets( ibuf );
}

