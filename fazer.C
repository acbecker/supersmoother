#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
//#include <fcntl.h>
//#include <sys/stat.h>
//#include <errno.h>
//#include <unistd.h>
#include <string.h>
#include <Macho.h>
#include <libsite.H>
#include <libclif.H>
#include "phase.H"
#include "fit/lcfit.H"
#include <curses.h>

/*---------------------------------------------------------------------------*
  fazer is designed to work with clc cached lightcurves.
 *---------------------------------------------------------------------------*/

enum fittype { None = 0, Linear, Cosine, Quadratic };
int main( int argc, char *argv[] ) {

   clifpar.load( );
   SodophotSet ss;
   CommonLc *rLc, *bLc, *tmpLc, lineLc1, lineLc2;
   LcPlot dplot, pplot, dpplot;
   DPhase rdp( clifpar ), bdp( clifpar );
   LcSet dlcs, plcs, dplcs, tlcs;
   ReimannPhaser rp( clifpar );
   CosineFit cosine( clifpar );
   LcMoment mom;
   
   char lcdir[512], path[512], inrec[512], outfi[512];
   char starid[64];
   char dev[256];
   char comment[256];
   char ylabel[256]; 
   double period, rper, bper;
   double hisig = 10.0, losig = 10.0, sigerr;
   double S, Sx, Sy, Sxx, S3x, S4x, Sxy, Sxxy;
   double tt, vv, ll;
   double ie2, lin, con, tunits;
   double minv, minp, dbuf;
   double binamp, binper, binphase0;
   float tlo, thi;
   float dtlo, dthi, dvlo, dvhi;
   float ptlo, pthi, pvlo, pvhi;
   double d, dd;
   double **A, **B;
   int field, tile, seq, chunk;
   int i, s, n, m, np, p, overlay;
   int c;
   fittype currentfit = None;
   int resflag = 0;

   
   FILE *fpin = 0;
   FILE *fpout = 0;
   FILE *fpid = 0;
   int single = 0;

   matrix( A, 3, 3 );
   matrix( B, 3, 1 );
   
   if( argc < 2 ) {
      fprintf( stderr, "Usage:\n\t%s field.tile.seq period\n",argv[0] );
      fprintf( stderr, "\t%s fts_period_file\n", argv[0] );
      return 1;
   } else if( argc == 2 ) {
      if( (fpin = fopen( argv[1], "r" ) ) == 0 ) {
	 fprintf( stderr, "Cannot open %s to read\n", argv[1] );
	 return 1;
      }
   } else {	
      sprintf( starid, "%s", argv[1] );
      if( sscanf( argv[2], "%lf", &period ) != 1 ) {
	 fprintf( stderr, "Cannot read period: %s\n", argv[2] );
	 return 1;
      }
      single = 1;
   }

   sprintf( outfi, "%s.notes", argv[0] );
   if( (fpid = fopen( outfi, "a" ) ) == 0  )
      fprintf(stderr,  "Cannot append to %s, continuing anyway\n ", outfi );
       
   if( getenv( "LCDIR" ) )
      sprintf( lcdir, "%s", getenv("LCDIR") );
   else	
      sprintf( lcdir, "%s/lcdir", getenv("HOME") );
   SodophotCuts sodcuts( clifpar, lcdir );

   
   /*--- curses stuff ---*/
   WINDOW *w = initscr( );
   cbreak( );
   //noecho( );

   /*--- make straight lines - could be prettier ---*/
   tlcs.append( tlcs.newlc( ) );
   tlcs.append( tlcs.newlc( ) );
   tlcs[0]->resize( 200 );
   tlcs[1]->resize( 200 );
   tlcs[0]->phtpkg = tlcs[1]->phtpkg = CommonLc::Theory;
   tlcs[0]->units = tlcs[1]->units = CommonLc::Linear;
   for( i = 0; i < 200; ++i ) {
      tlcs[0]->date[i] = tlcs[1]->date[i] = 10.0 * i;
      tlcs[0]->value[i] = tlcs[1]->value[i] = 0.0;
      tlcs[0]->error[i] = tlcs[1]->error[i] = 0.0;
   }

   while( single || fgets( inrec, 512, fpin ) ) {	// loop over stars
      if( fpin ) {
	 if( *inrec == '#' || *inrec == '!' ||
	     sscanf( inrec, "%s %lf", starid, &period ) != 2 ) 
	    continue;
      } else
	 single = 0;		// only works the first time through
      
      /*--- recycle old data ---*/
      dlcs.recycle( );

      /*--- plcs and dplcs lightcurves are handled individually ---*/
      while( plcs[0] )
	 plcs.unhook( plcs[0] );
      while( dplcs[0] )
	 dplcs.unhook( dplcs[0] );

      /*--- load and merge new data ---*/
      lcdir_sodlc_file( lcdir, starid, 'p', path );
      if( !(ss.open( path ) || sodcuts.clean( dlcs, ss, 0 )) ) {
	 chunk = ss.star[0].wcid[Red];		// success
      }
      lcdir_sodlc_file( lcdir, starid, 'o', path );
      if( !(ss.open( path ) || sodcuts.clean( dlcs, ss, 0 )) ) {
	 chunk = ss.star[0].wcid[Red];		// success
      }
      LcMerge( dlcs );
      if( dlcs.count( ) == 0  ) {
	 goto nextstar;				// failure
      }
      for( i = 0; i < dlcs.count( ); ++i ) {
	 dlcs[i]->convert( CommonLc::Linear );
	 mom.eval( *dlcs[i] );
	 dlcs[i]->renorm( mom.average( ), 0.0 );
      }
      
      
      /*--- find red and blue data lightcurves ---*/
      rLc = bLc = 0;
      if( dlcs[0] ) {
	 if( Passband::info((*dlcs[0]).filter)->id == Passband::R_MACHO ) 
	    rLc = dlcs[0];
	 else if(Passband::info((*dlcs[0]).filter)->id == Passband::V_MACHO )
	    bLc = dlcs[0];
      }
      if( dlcs[1] ) {
	 if( Passband::info((*dlcs[1]).filter)->id == Passband::R_MACHO ) 
	    rLc = dlcs[1];
	 else if( Passband::info((*dlcs[1]).filter)->id == Passband::V_MACHO )
	    bLc = dlcs[1];
      }
      if( (bLc ? bLc->nobs : 0) + (rLc ? rLc->nobs : 0) < 50 ) {
	 printw("starid: %s period: %13.8f   Nred: %4d   Nblue: %4d > ",
		starid,  period, rLc ? rLc->nobs : 0, bLc ? bLc->nobs : 0 );
	 printw("Continuing - not enough observations\n");
	 continue;
      }

      /*--- plot data ---*/
      dplot.config( dlcs );
      sprintf( dplot.xlabel, "JD - %.5lf (P = %.8lf d)",
	       dlcs[0]->refJD, period );
      dplot.draw( dlcs );

      /*--- calculate ss curve and fit dphase curve ---*/
      if( rLc && rLc->nobs > 30 && !rdp.derive( *rLc, period ) ) {
	 plcs.append( &rdp.pLc );
	 plcs.append( &rdp.ssLc );
	 dplcs.append( &rdp.dpLc );
	 lc_hdr_cpy( *rLc, *tlcs[0] );
	 tlcs[0]->phtpkg = CommonLc::Theory;
	 tlcs[0]->units = CommonLc::Linear;
	 dplcs.append( tlcs[0] );
      } 
      if( bLc && bLc->nobs > 30 && !bdp.derive( *bLc, period ) ) {
	 plcs.append( &bdp.pLc );
	 plcs.append( &bdp.ssLc );
	 dplcs.append( &bdp.dpLc );
	 lc_hdr_cpy( *bLc, *tlcs[1] );
	 tlcs[1]->phtpkg = CommonLc::Theory;
	 tlcs[1]->units = CommonLc::Linear;
	 dplcs.append( tlcs[1] );
      }
      /*--- set phase plot defaults ---*/  
      pplot.config( plcs );
      sprintf( pplot.xlabel, "Phase/2\\pi\\ \\ (P = %.8lf d)", period );
      pplot.xlimlo = 0.0;
      pplot.xlimhi = 1.0;

      /*--- set phase plot defaults ---*/
      dpplot.config( dplcs );
      if( (tunits = rdp.ok ? rdp.tunits : bdp.tunits) == 1.0 ) 
	 sprintf( ylabel, "Time Delay (\\Delta t/P)" );
      else
	 sprintf( ylabel, "Time Delay (AU/c)" );
      sprintf( dpplot.xlabel, "JD - %.5lf (P = %.8lf d, B = %d)",
	       dlcs[0]->refJD, period, bdp.obsperfit );
      sprintf( dpplot.ylabel, "%s", ylabel );
      dpplot.xlimlo = dplot.xlimlo;	// keep same limits as data
      dpplot.xlimhi = dplot.xlimhi;	// keep same limits as data

      while( 1 ) {	/*--- print menu and enter options loop ---*/
	 
	 clear( );
	 printw("starid: %-17s   period: %.8f d\n", starid, period );
	 if( rdp.ok ) printw("Red   N: %-4d  fom: %8.4f  chi2: %8.4f\n",
			     rLc->nobs, rdp.fom( ), rdp.chi2( ) ); 
	 if( bdp.ok ) printw("Blue  N: %-4d  fom: %8.4f  chi2: %8.4f\n",
			     bLc->nobs, bdp.fom( ), bdp.chi2( ) );
	 printw("\n");
	 printw("D  Data lightcurve\n");
	 printw("P  Phased lightcurve\n");
	 printw("O  Time delay lightcurve\n");
	 //printw("T  Time limits\n" );
	 printw("M  Data minus smoothed fit (residual toggle)\n" );
	 printw("N  Next lightcurve\n" );
	 printw("S  Save starID, period and comment\n" );
	 printw("Q  Quit\n" );
	 printw("B  Bin size\n");
	 printw("C  Clean data (sigma clip phased curve)\n" );
	 printw("E  Print or write eps files (device = file)\n");
	 printw("F  Fit sin to time delay curve\n" ); 
	 printw("G  Refresh data\n" ); 
	 printw("J  Redefine phase 0 as minimum light (JD offset)\n"); 
	 printw("L  Refine period estimate (linear time delay fit)\n" );
	 printw("R  Recalculate period (Reimann phaser)\n");
	 //printw("U  Toggle time delay units\n");
	 //printw("V  Toggle full data overlay\n");
	 printw("W  Write time delay data file\n" );
	 //printw("Z  Insert fake signal\n");
	 printw("\n> " );
	 c = getch( );

	 switch( toupper( c ) ) {

	 case 'D' :	// plot raw data
	    dplot.draw( dlcs );
	    break;
	    
	 case 'P' :	// plot phased data w/ ss curves
	    pplot.draw( plcs );
	    break;
	    
	 case 'O' :	// plot time delay date
	    dpplot.draw( dplcs );
	    break;
	    
	 case 'V' :	// toggle ss overlay for full timeseries
	    if( overlay ) {
	       n = 0;
	       while( n < dlcs.count( ) ) {
		  if( dlcs[n]->phtpkg == CommonLc::Theory )
		     dlcs.recycle( dlcs[n] );
		  else
		     ++n;
	       }
	       overlay = 0;
	       dplot.draw( dlcs );
	       break;
	    }
	    dd = period / 100;
	    n = int( (dthi - dtlo) / dd ) + 1;
	    if( rdp.ok ) {
	       tmpLc = dlcs.newlc( );
	       lc_hdr_cpy( *rLc, *tmpLc );
	       tmpLc->phtpkg = CommonLc::Theory;
	       tmpLc->resize( n );
	       for( m = 0; m < n; ++m ) {
		  tmpLc->date[m] = dtlo + m * dd;
		  tmpLc->value[m] = rdp.ss( dtlo + m * dd );
		  tmpLc->error[m] = 0;
	       }
	       dlcs.append( tmpLc );
	       overlay = 1;
	    }
	    if( bdp.ok ) {
	       tmpLc = dlcs.newlc( );
	       lc_hdr_cpy( *bLc, *tmpLc );
	       tmpLc->phtpkg = CommonLc::Theory;
	       tmpLc->resize( n );
	       for( m = 0; m < n; ++m ) {
		  tmpLc->date[m] = dtlo + m * dd;
		  tmpLc->value[m] = bdp.ss( dtlo + m * dd );
		  tmpLc->error[m] = 0;
	       }
	       dlcs.append( tmpLc );
	       overlay = 1;
	    }
	    dplot.draw( dlcs );
	    break;

	 case 'S' :
	    fprintf( fpid, "%-17s  %16.8f %16.6f",
		     starid, period, dlcs[0]->refJD );
	    printw("\nComment (no spaces) > ");
	    scanw( "%s", comment );
	    if( *comment )
	       fprintf( fpid, "\t#%s\n", comment );
	    else
	       fprintf( fpid, "\n" );
	    fflush( fpid );
	    break;

	 case 'M' :
	    if( !resflag ) {
	       if( rLc ) for( i = 0; i < rLc->nobs; ++i )
		  rLc->value[i] -= rdp.ss( rLc->date[i] );
	       if( bLc ) for( i = 0; i < bLc->nobs; ++i )
		  bLc->value[i] -= bdp.ss( bLc->date[i] );
	       dplot.config( dlcs );
	       dplot.draw( dlcs );
	       resflag = 1;
	    } else {
	       if( rLc ) for( i = 0; i < rLc->nobs; ++i )
		  rLc->value[i] += rdp.ss( rLc->date[i] );
	       if( bLc ) for( i = 0; i < bLc->nobs; ++i )
		  bLc->value[i] += bdp.ss( bLc->date[i] );
	       dplot.config( dlcs );
	       dplot.draw( dlcs );
	       resflag = 0;
	    }
	    break;
	    
	 case 'N' :
	    currentfit = None;
	    goto nextstar;

	 case 'Q' :
	    printw("Quit y/n > ");
	    c = getch( );
	    if( toupper( c ) != 'Y' )
	       break;
	    endwin( );

	    /*--- plcs and dplcs lightcurves are handled by DPhase class ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );

	    if( fpid )
	       fclose( fpid );
	    return 0;

	 case 'T' :
	    printw("\nLow date > " );
	    scanw("%lf", &dplot.xlimlo );
	    printw("High date > " );
	    scanw("%f", &dplot.xlimhi );
	    dpplot.xlimlo = dplot.xlimlo;
	    dpplot.xlimhi = dplot.xlimhi;
	    dplot.draw( dlcs );
	    break;

	 case 'B' :	// redefine obsperfit in oc curve
	    printw("\nObservations per fit (%d,%d) > ",
		   rdp.obsperfit, bdp.obsperfit );
	    scanw("%d", &i );
	    if( i < 0 ||i > (rdp.ndata + bdp.ndata) / 2 ) {
	       printw("Invalid bin size (%d)\n");
	       getch( );
	       continue;
	    }
	    rdp.obsperfit = bdp.obsperfit = i;
	    sprintf( dpplot.xlabel, "JD - %.5lf (P = %.8lf d, B = %d)",
		     dlcs[0]->refJD, period, bdp.obsperfit );

	    /*--- fold, calculate ss curve and fit dphase curve ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );
	    if( rLc && rLc->nobs > 30 && !rdp.derive( *rLc, period ) ) {
	       plcs.append( &rdp.pLc );
	       plcs.append( &rdp.ssLc );
	       dplcs.append( &rdp.dpLc );
	       dplcs.append( tlcs[0] );
	    } 
	    if( bLc && bLc->nobs > 30 && !bdp.derive( *bLc, period ) ) {
	       plcs.append( &bdp.pLc );
	       plcs.append( &bdp.ssLc );
	       dplcs.append( &bdp.dpLc );
	       dplcs.append( tlcs[1] );
	    }
	    LcDataLimits( dplcs, tlo, thi, pvlo, pvhi );
	    pvhi += 0.2 * (pvhi - pvlo );
	    dpplot.ylimlo = pvlo;
	    dpplot.ylimhi = pvhi;
	    sprintf( dpplot.ylabel, "%s", ylabel );
	    dpplot.draw( dplcs );
	    break;
	    
	 case 'C' :	// clean by sigma clipping from ss curve 
	    printw("\nHigh sigma limit (%3.1f) > ", losig );	// dimnitudes
	    c = getch( );
	    if( isdigit( c ) )
	       losig = c - 48;
	    printw("\nLow sigma limit(%3.1f)  > ", hisig );	// dimnitudes
	    c = getch( );
	    if( isdigit( c ) )
	       hisig = c - 48;
	    if( rLc ) {
	       tmpLc = dlcs.newlc( );
	       lc_hdr_cpy( *rLc, *tmpLc );
	       tmpLc->resize( rLc->nobs );
	       for( m = n = 0; n < rLc->nobs; ++n ) {
		  sigerr = (rLc->value[n] - rdp.ssLc.value[rdp.dindex[n]])/
		     rLc->error[n];
		  if( sigerr > -losig && sigerr < hisig ) {
		     tmpLc->date[m] = rLc->date[n];
		     tmpLc->value[m] = rLc->value[n];
		     tmpLc->error[m] = rLc->error[n];
		     tmpLc->obsid[m] = rLc->obsid[n];
		     ++m;
		  }
	       }
	       dlcs.recycle( rLc );
	       rLc = tmpLc;
	       rLc->nobs = m;
	       dlcs.append( rLc );
	    }

	    if( bLc ) {
	       tmpLc = dlcs.newlc( );
	       lc_hdr_cpy( *bLc, *tmpLc );
	       tmpLc->resize( bLc->nobs );
	       for( m = n = 0; n < bLc->nobs; ++n ) {
		  sigerr = (bLc->value[n] - bdp.ssLc.value[bdp.dindex[n]])/
		     bLc->error[n];
		  if( sigerr > -losig && sigerr < hisig ) {
		     tmpLc->date[m] = bLc->date[n];
		     tmpLc->value[m] = bLc->value[n];
		     tmpLc->error[m] = bLc->error[n];
		     tmpLc->obsid[m] = bLc->obsid[n];
		     ++m;
		  }
	       }
	       dlcs.recycle( bLc );
	       bLc = tmpLc;
	       bLc->nobs = m;
	       dlcs.append( bLc );
	    }
	    /*--- plot data ---*/
	    LcDataLimits( dlcs, dtlo, dthi, dvlo, dvhi );
	    dvlo -= 0.2 * ( dvhi - dvlo );		// for magnitudes

	    /*--- fold, calculate ss curve and fit dphase curve ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );
	    if( rLc && rLc->nobs > 30 && !rdp.derive( *rLc, period ) ) {
	       plcs.append( &rdp.pLc );
	       plcs.append( &rdp.ssLc );
	       dplcs.append( &rdp.dpLc );
	       dplcs.append( tlcs[0] );
	    } 
	    if( bLc && bLc->nobs > 30 && !bdp.derive( *bLc, period ) ) {
	       plcs.append( &bdp.pLc );
	       plcs.append( &bdp.ssLc );
	       dplcs.append( &bdp.dpLc );
	       dplcs.append( tlcs[1] );
	    }
	    LcDataLimits( dplcs, tlo, thi, pvlo, pvhi );
	    pvhi += 0.2 * (pvhi - pvlo );
	    dpplot.ylimlo = pvlo;
	    dpplot.ylimhi = pvhi;

	    pplot.draw( plcs );
	    break;

	 case 'E' :	// write eps files / print
	    printw("\ndevice > ");
	    scanw("%s", dev );
	    printw("\n");
	    if( strncmp( dev, "file", 4 ) == 0 ) {
	       sprintf( outfi, "%s.dat.eps", starid );
	       dplot.print( dlcs, "postencap",  outfi );
	       printw("Writing %s\n", outfi );
	       sprintf( outfi, "%s.per.eps", starid );
	       pplot.print( plcs, "postencap",  outfi );
	       printw("Writing %s\n", outfi );
	       sprintf( outfi, "%s.oc.eps", starid );
	       dpplot.print( dplcs, "postencap",  outfi );
	       printw("Writing %s\n", outfi );
	    } else {
	       dplot.print( dlcs, dev  );
	       pplot.print( plcs, dev );
	       dpplot.print( dplcs, dev );
	    }
	    printw("> " );
	    getch( );
	    break;

	 case 'G' :

	    /*--- recycle old data ---*/
	    dlcs.recycle( );

	    /*--- plcs and dplcs lightcurves are handled individually ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );

	    /*--- load and merge new data ---*/
	    lcdir_sodlc_file( lcdir, starid, 'p', path );
	    if( !(ss.open( path ) || sodcuts.clean( dlcs, ss, 0 )) ) {
	       chunk = ss.star[0].wcid[Red];		// success
	    }
	    lcdir_sodlc_file( lcdir, starid, 'o', path );
	    if( !(ss.open( path ) || sodcuts.clean( dlcs, ss, 0 )) ) {
	       chunk = ss.star[0].wcid[Red];		// success
	    }
	    LcMerge( dlcs );
	    LcDataLimits( dlcs, dtlo, dthi, dvlo, dvhi );
	    dvlo -= 0.2 * ( dvhi - dvlo );		// for magnitudes

	    /*--- fold, calculate ss curve and fit dphase curve ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );
	    if( rLc && rLc->nobs > 30 && !rdp.derive( *rLc, period ) ) {
	       plcs.append( &rdp.pLc );
	       plcs.append( &rdp.ssLc );
	       dplcs.append( &rdp.dpLc );
	       dplcs.append( tlcs[0] );
	    } 
	    if( bLc && bLc->nobs > 30 && !bdp.derive( *bLc, period ) ) {
	       plcs.append( &bdp.pLc );
	       plcs.append( &bdp.ssLc );
	       dplcs.append( &bdp.dpLc );
	       dplcs.append( tlcs[1] );
	    }

	    /*--- resize dpplot as values may have changed ---*/ 
	    LcDataLimits( dplcs, tlo, thi, pvlo, pvhi );
	    pvhi += 0.2 * (pvhi - pvlo );
	    dpplot.ylimlo = pvlo;
	    dpplot.ylimhi = pvhi;

	    dplot.draw( dlcs );
	    break;
	    
	 case 'J' :
	    if( !plcs.count( ) )
	       break;
	    if( resflag ) {
	       if( rLc ) for( i = 0; i < rLc->nobs; ++i )
		  rLc->value[i] += rdp.ss( rLc->date[i] );
	       if( bLc ) for( i = 0; i < bLc->nobs; ++i )
		  bLc->value[i] += bdp.ss( bLc->date[i] );
	       dplot.config( dlcs );
	       resflag = 0;
	    }
	    minv = -Infinity;		// dimnitudes of course
	    minp = 0;
	    for( m = 0; m < plcs.count( ); ++m ) {
	       if( plcs[m]->phtpkg != CommonLc::Theory )	// us ss curves
		  continue;
	       for( i = 0; i < plcs[m]->nobs; ++i ) {
		  if( plcs[m]->value[i] > minv ) {
		     minv = plcs[m]->value[i];
		     minp = plcs[m]->date[i];
		  }
	       }
	    }
	    if( isinf( minv ) )
	       break;
	    minp *= period;
	    for( m = 0; m < dlcs.count( ); ++m ) {
	       dlcs[m]->refJD += minp;
	       for( i = 0; i < dlcs[m]->nobs; ++i ) 
		  dlcs[m]->date[i] -= minp;
	    }
	    sprintf( dplot.xlabel, "JD - %.5lf (P = %.8lf d)",
		     dlcs[0]->refJD, period );
	    sprintf( dpplot.xlabel, "JD - %.5lf (P = %.8lf d, B = %d)",
		     dlcs[0]->refJD, period, rdp.obsperfit );
	    
	    /*--- fold, calculate ss curve and fit dphase curve ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );
	    if( rLc && rLc->nobs > 30 && !rdp.derive( *rLc, period ) ) {
	       plcs.append( &rdp.pLc );
	       plcs.append( &rdp.ssLc );
	       dplcs.append( &rdp.dpLc );
	       dplcs.append( tlcs[0] );
	    } 
	    if( bLc && bLc->nobs > 30 && !bdp.derive( *bLc, period ) ) {
	       plcs.append( &bdp.pLc );
	       plcs.append( &bdp.ssLc );
	       dplcs.append( &bdp.dpLc );
	       dplcs.append( tlcs[1] );
	    }
	    LcDataLimits( dplcs, tlo, thi, pvlo, pvhi );
	    pvhi += 0.2 * (pvhi - pvlo );
	    dpplot.ylimlo = pvlo;
	    dpplot.ylimhi = pvhi;

	    pplot.draw( plcs );
	    break;
	    
	 case 'L' :	// improve period estimate with linear fit to dplc
	    S = Sx = Sy = Sxx = Sxy = 0;
	    for( p = m = 0; m < dplcs.count( ); ++m ) {
	       if( dplcs[m]->phtpkg == CommonLc::Theory )
		  continue;
	       for( n = 0; n < dplcs[m]->nobs; ++n ) {
		  ++p;
		  S += ie2 = pow(dplcs[m]->error[n], -2.0 );
		  Sx += dplcs[m]->date[n] * ie2;
		  Sy += dplcs[m]->value[n] * ie2;
		  Sxx += dplcs[m]->date[n] * dplcs[m]->date[n] * ie2;
		  Sxy += dplcs[m]->date[n] * dplcs[m]->value[n] * ie2;
	       }
	    }
	    sprintf( dpplot.plabel, "star ID: %s", starid );
	    tunits = rdp.ok ? rdp.tunits : bdp.tunits;
	    lin = (S * Sxy - Sx * Sy ) / (S * Sxx - Sx * Sx );
	    con = (Sxx * Sy - Sx * Sxy ) / (S * Sxx - Sx * Sx );
	    for( i = 0; i < 200; ++i ) {
	       tlcs[0]->value[i] = tlcs[1]->value[i] =
		  con + lin * tlcs[0]->date[i];
	    }
	    currentfit = Linear;
	    printw( "\nfit = %.6g + t * %.6g for %d points\n",
		    lin, con, p );
	    //lin *= -pow(period, 2.0) / tunits;
	    lin *= -period / tunits;
	    sprintf( dpplot.ylabel, "%s", ylabel );
	    dpplot.draw( dplcs );
	    printw( "Correct period by %.8lf d y/n > ", lin );
	    c = getch( );
	    printw("\n");

	    for( i = 0; i < 200; ++i )	// reset lines to zero
	       tlcs[0]->value[i] = tlcs[1]->value[i] = 0;
	    if( toupper(c) != 'Y' )
	       break;

	    period += lin;
	    sprintf( dplot.xlabel, "JD - %.5lf (P = %.8lf d)",
		     dlcs[0]->refJD, period );
	    sprintf( pplot.xlabel, "Phase/2\\pi\\ \\ (P = %.8lf d)", period );
	    sprintf( dpplot.xlabel, "JD - %.5lf (P = %.8lf d, B = %d)",
		     dlcs[0]->refJD, period, rdp.obsperfit );

	    /*--- fold, calculate ss curve and fit dphase curve ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );
	    if( rLc && rLc->nobs > 30 && !rdp.derive( *rLc, period ) ) {
	       plcs.append( &rdp.pLc );
	       plcs.append( &rdp.ssLc );
	       dplcs.append( &rdp.dpLc );
	       dplcs.append( tlcs[0] );
	    } 
	    if( bLc && bLc->nobs > 30 && !bdp.derive( *bLc, period ) ) {
	       plcs.append( &bdp.pLc );
	       plcs.append( &bdp.ssLc );
	       dplcs.append( &bdp.dpLc );
	       dplcs.append( tlcs[1] );
	    }
	    LcDataLimits( dplcs, tlo, thi, pvlo, pvhi );
	    pvhi += 0.2 * (pvhi - pvlo );
	    dpplot.ylimlo = pvlo;
	    dpplot.ylimhi = pvhi;

	    sprintf( dpplot.ylabel, "%s", ylabel );
	    dpplot.draw( dplcs );
	    break;
	    
	 case 'F' :
	    dplcs.unhook( tlcs[0] );
	    dplcs.unhook( tlcs[1] );
	    currentfit = None;
	    printw(" (C)osine  (Q)uadratic > ");
	    c = getch( );
	    switch( toupper( c ) ) {
	    case 'C' :
	       cosine.init( dplcs );
	       printw( "\nAmplitude > " );
	       scanw( "%lf", &dbuf );
	       cosine.init_par( 0, dbuf, GeneralLcFit::Fit );
	       printw( "Period (days) > " );
	       scanw( "%lf", &dbuf );
	       cosine.init_par( 1, dbuf, GeneralLcFit::Fit );
	       printw( "Zero Phase (days) > " );
	       scanw( "%lf", &dbuf );
	       cosine.init_par( 2, dbuf, GeneralLcFit::Fit );
	       if( cosine.ifit(0, 0, 0 ) == SiteErr::errNone ) {
		  sprintf( dpplot.plabel,
			   "star ID: %s orbit:  A=%.2f P=%.1f \\theta=%.1f",
			   starid,cosine.par[0], cosine.par[1], cosine.par[2]);
		  cosine.eval(tlcs[0]->date,tlcs[0]->value,0,tlcs[0]->nobs );
		  cosine.eval(tlcs[1]->date,tlcs[1]->value,0,tlcs[1]->nobs );
		  dplcs.append( tlcs[0] );
		  dplcs.append( tlcs[1] );
		  dpplot.draw( dplcs );
		  printw("Log fit (y/n) >" );
		  c = getch( );
		  if( toupper( c ) == 'Y' )
		     fprintf(fpid,"%-17s %12.8f "
			     "1 %3d %.5g %.3g %.5g %.3g %.5g %.3g %.3g\n",
			     starid, period, bdp.obsperfit, 
			     cosine.par[0], sqrt(cosine.covar[0][0]),
			     cosine.par[1], sqrt(cosine.covar[1][1]),
			     cosine.par[2], sqrt(cosine.covar[2][2]),
			     cosine.Chi2 / cosine.DOF );
		  currentfit = Cosine;
	       } else {
		  printw( "Could not fit > ");
		  getch( );
	       }
	       break;
	    case 'Q' :
	       S = Sx = Sxx = S3x = S4x = Sy = Sxy = Sxxy = 0;
	       for( p = m = 0; m < dplcs.count( ); ++m ) {
		  if( dplcs[m]->phtpkg == CommonLc::Theory )
		     continue;
		  for( n = 0; n < dplcs[m]->nobs; ++n ) {
		     ++p;
		     S += ie2 = pow(dplcs[m]->error[n], -2.0 );
		     Sx += ll = ((tt = dplcs[m]->date[n]) * ie2);
		     Sy += (vv = dplcs[m]->value[n]) * ie2;
		     Sxy += ll * vv;
		     Sxx += ll *= tt;
		     Sxxy += ll * vv;
		     S3x += ll *= tt;
		     S4x += ll *= tt;
		  }
	       }
	       A[0][0] = S;
	       A[1][0] = A[0][1] = Sx;
	       A[1][1] = A[2][0] = A[0][2] = Sxx;
	       A[2][1] = A[1][2] = S3x;
	       A[2][2] = S4x;
	       B[0][0] = Sy;
	       B[1][0] = Sxy;
	       B[2][0] = Sxxy;
	       cgaussj( A, 3, B, 1 );
	       if( 1 ) {
		  printw("\n%.3g + %.4g t + %.5g t^2   ",
			 B[0][0], B[1][0], B[2][0] );	
		  printw("dP/dt = %.5lg\n", -2.0 * period * B[2][0] / tunits );
		  sprintf( dpplot.plabel,  "star ID: %s\\ \\ \\ dP/dt = %.5g",
			   starid, -2.0 * period * B[2][0] / tunits );
		  for( n = 0; n < tlcs[0]->nobs; ++n )
		     tlcs[0]->value[n] = B[0][0] + B[1][0]* tlcs[0]->date[n] +
			B[2][0] * tlcs[0]->date[n] * tlcs[0]->date[n];
		  for( n = 0; n < tlcs[1]->nobs; ++n )
		     tlcs[1]->value[n] = B[0][0] + B[1][0]* tlcs[1]->date[n] +
			B[2][0] * tlcs[1]->date[n] * tlcs[1]->date[n];
		  dplcs.append( tlcs[0] );
		  dplcs.append( tlcs[1] );
		  dpplot.draw( dplcs );
		  printw("Log fit (y/n) >" );
		  c = getch( );
		  if( toupper( c ) == 'Y' )
		     fprintf(fpid,"%-17s %12.8f 2 %3d %.4g %.4e %.4e %.4e\n",
			     starid, period, bdp.obsperfit, 
			     B[0][0], B[1][0], B[2][0],
			     -2.0 * period * B[2][0] / tunits );
		  currentfit = Quadratic;
	       } else {
		  printw( "Could not fit > ");
		  getch( );
	       }
	       break;
	    default :
	       printw("Invalid fit type > " );
	       getch( ); 
	       break;
	    }
	    break;
	    
	 case 'R' :
	    printw("\nReimannPhaser.MINP > " );
	    scanw( "%lf", &rp.MINP );
	    printw("ReimannPhaser.MAXP > " );
	    scanw( "%lf", &rp.MAXP );
	    printw("\n" );
	    np = 0;
	    rper = bper = period = 0;
	    if( rLc ) {
	       if( rp.phase( *rLc ) ) {
		  printw("Could not determine red period\n");
	       } else {
		  ++np;
		  period = rper = 1.0/rp.FF[0];
	       }
	    }
	    if( bLc ) {
	       if( rp.phase( *bLc ) ) {
		  printw("Could not determine blue period\n");
	       } else {
		  ++np;
		  period += bper = 1.0/rp.FF[0];
	       }
	    }
	    period /= np;
	    printw("Red period: %.8lf  Blue period: %.8lf\n",
		   rper ? rper : NaN, bper ? bper : NaN );
	    printw("Use average %.8f  y/n >", period );
	    c = getch( );
	    if( toupper( c ) != 'Y' ) {
	       printw("Enter period > " );
	       scanw( "%lf", &period );
	       printw("\n");
	    }

	    sprintf( dplot.xlabel, "JD - %.5lf (P = %.8lf d)",
		     dlcs[0]->refJD, period );
	    sprintf( pplot.xlabel, "Phase/2\\pi\\ \\ (P = %.8lf d)", period );
	    sprintf( dpplot.xlabel, "JD - %.5lf (P = %.8lf d, B = %d)",
		     dlcs[0]->refJD, period, rdp.obsperfit );

	    /*--- fold, calculate ss curve and fit dphase curve ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );
	    if( rLc && rLc->nobs > 30 && !rdp.derive( *rLc, period ) ) {
	       plcs.append( &rdp.pLc );
	       plcs.append( &rdp.ssLc );
	       dplcs.append( &rdp.dpLc );
	       dplcs.append( tlcs[0] );
	    } 
	    if( bLc && bLc->nobs > 30 && !bdp.derive( *bLc, period ) ) {
	       plcs.append( &bdp.pLc );
	       plcs.append( &bdp.ssLc );
	       dplcs.append( &bdp.dpLc );
	       dplcs.append( tlcs[1] );
	    }
	    LcDataLimits( dplcs, tlo, thi, pvlo, pvhi );
	    pvhi += 0.2 * (pvhi - pvlo );
	    dpplot.ylimlo = pvlo;
	    dpplot.ylimhi = pvhi;
	    pplot.draw( plcs );
	    break;

	 case 'Z' :	// insert test signal
	    double tspan;
	    printw("\nBinary period (days) > " );
	    scanw( "%lf", &binper );
	    printw("Amplitude (AU) > " );
	    scanw( "%lf ", &binamp );
	    printw("Time of closest approach (days) > " );
	    scanw( "%lf ", &binphase0 );
	    binamp /= 173.264498255;		// days
	    
	    if( rLc ) {
	       printw("Defacing red curve\n");
	       for( i = 0; i < rLc->nobs; ++i ) 
		  rLc->date[i] -= binamp *
		     cos(2.0*M_PI*(rLc->date[i]-binphase0 )/binper);
	    }	
	    if( bLc ) {
	       printw("Defacing blue curve\n");
	       for( i = 0; i < bLc->nobs; ++i ) 
		  bLc->date[i] -= binamp *
		     cos(2.0*M_PI*(bLc->date[i]-binphase0 )/binper);
	    }
	    getch( );

	    /*--- fold, calculate ss curve and fit dphase curve ---*/
	    while( plcs[0] )
	       plcs.unhook( plcs[0] );
	    while( dplcs[0] )
	       dplcs.unhook( dplcs[0] );
	    if( rLc && rLc->nobs > 30 && !rdp.derive( *rLc, period ) ) {
	       plcs.append( &rdp.pLc );
	       plcs.append( &rdp.ssLc );
	       dplcs.append( &rdp.dpLc );
	       dplcs.append( tlcs[0] );
	    } 
	    if( bLc && bLc->nobs > 30 && !bdp.derive( *bLc, period ) ) {
	       plcs.append( &bdp.pLc );
	       plcs.append( &bdp.ssLc );
	       dplcs.append( &bdp.dpLc );
	       dplcs.append( tlcs[1] );
	    }
	    LcDataLimits( dplcs, tlo, thi, pvlo, pvhi );
	    pvhi += 0.2 * (pvhi - pvlo );
	    dpplot.ylimlo = pvlo;
	    dpplot.ylimhi = pvhi;
	    sprintf( dpplot.ylabel, "%s", ylabel );
	    dpplot.draw( dplcs );
	    break;

	 case 'U' :
	    if( tunits == 1.0 ) {
	       sprintf( ylabel, "Time Delay (AU/c)" );
	       rdp.tunits = bdp.tunits = tunits =
		  86400*2.997924e10/1.495979e13;	// days --> AU
	       for( m = 0; m < dplcs.count( ); ++m ) 
		  dplcs[m]->renorm( 1.0/tunits, 0.0 );
	       dpplot.ylimlo *= tunits;
	       dpplot.ylimhi *= tunits;
	    } else {
	       sprintf( ylabel, "Time Delay (days)" );
	       for( m = 0; m < dplcs.count( ); ++m )
		  dplcs[m]->renorm( tunits, 0.0 );
	       dpplot.ylimlo /= tunits;
	       dpplot.ylimhi /= tunits;
	       rdp.tunits = bdp.tunits = tunits = 1.0;
	    }

	    sprintf( dpplot.ylabel, "%s", ylabel );
	    dpplot.draw( dplcs );
	    break;

	 case 'W' :		// write time delay data file
	    if( rLc && rdp.ok ) {
	       sprintf( outfi, "%s.rdt", starid );
	       if( (fpout = fopen( outfi, "w" ) ) == 0 ) {
		  printw("\nCannot open %s to write\n", starid );
	       } else {
		  switch( currentfit ) {
		  case Linear :
		     fprintf( fpout, "# P = %.7g, bin = %d\n",
			      period, rdp.obsperfit );
		     fprintf( fpout, "# %.6g + T * %.6g\n", lin, con );
		     break;
		  case Cosine :
		     fprintf( fpout, "# P = %.7g, bin = %d\n",
			      period, rdp.obsperfit );
		     fprintf( fpout, "# %.4g cos( 2_PI * (T - %.4g) / %.4g)\n",
			      cosine.par[0], cosine.par[2], cosine.par[1]);
		     break;
		  case Quadratic :
		     fprintf( fpout, "# P = %.7g, bin = %d\n",
			      period, rdp.obsperfit );
		     fprintf(fpout, "# %.4g + %.4g T + %.5g T^2\n",
			    B[0][0], B[1][0], B[2][0] );	
		     fprintf( fpout, "# dP/dt = %.5lg\n",
			      -2.0 * period * B[2][0] / tunits );
		     break;
		  default :
		     break;
		  }
		  rdp.dpLc.print( fpout );
		  fclose( fpout );
	       }
	       sprintf( outfi, "%s.rlc", starid );
	       if( (fpout = fopen( outfi, "w" ) ) == 0 ) {
		  printw("\nCannot open %s to write\n", starid );
	       } else {
		  rLc->print( fpout ); 
		  fclose( fpout );
	       }
	       sprintf( outfi, "%s.rp", starid );
	       if( (fpout = fopen( outfi, "w" ) ) == 0 ) {
		  printw("\nCannot open %s to write\n", starid );
	       } else {
		  rdp.pLc.print( fpout ); 
		  fclose( fpout );
	       }
	       sprintf( outfi, "%s.rss", starid );
	       if( (fpout = fopen( outfi, "w" ) ) == 0 ) {
		  printw("\nCannot open %s to write\n", starid );
	       } else {
		  rdp.ssLc.print( fpout ); 
		  fclose( fpout );
	       }
	    }
	       
	    if( bLc && bdp.ok ) {
	       sprintf( outfi, "%s.bdt", starid );
	       if( (fpout = fopen( outfi, "w" ) ) == 0 ) {
		  printw("Cannot open %s to write\n", starid );
	       } else {
		  switch( currentfit ) {
		  case Linear :
		     fprintf( fpout, "# P = %.7g, bin = %d\n",
			      period, bdp.obsperfit );
		     fprintf( fpout, "# %.6g + T * %.6g\n", lin, con );
		     break;
		  case Cosine :
		     fprintf( fpout, "# P = %.7g, bin = %d\n",
			      period, bdp.obsperfit );
		     fprintf( fpout, "# %.4g cos( 2_PI * (T - %.4g) / %.4g)\n",
			      cosine.par[0], cosine.par[2], cosine.par[1]);
		     break;
		  case Quadratic :
		     fprintf( fpout, "# P = %.7g, bin = %d\n",
			      period, bdp.obsperfit );
		     fprintf(fpout, "# %.4g + %.4g T + %.5g T^2\n",
			    B[0][0], B[1][0], B[2][0] );	
		     fprintf( fpout, "# dP/dt = %.5lg\n",
			      -2.0 * period * B[2][0] / tunits );
		     break;
		  default :
		     break;
		  }
		  bdp.dpLc.print( fpout );
		  fclose( fpout );
	       }
	       sprintf( outfi, "%s.blc", starid );
	       if( (fpout = fopen( outfi, "w" ) ) == 0 ) {
		  printw("\nCannot open %s to write\n", starid );
		  endwin( );
	       } else {
		  bLc->print( fpout ); 
		  fclose( fpout );
	       }
	       sprintf( outfi, "%s.bp", starid );
	       if( (fpout = fopen( outfi, "w" ) ) == 0 ) {
		  printw("\nCannot open %s to write\n", starid );
	       } else {
		  bdp.pLc.print( fpout ); 
		  fclose( fpout );
	       }
	       sprintf( outfi, "%s.bss", starid );
	       if( (fpout = fopen( outfi, "w" ) ) == 0 ) {
		  printw("\nCannot open %s to write\n", starid );
	       } else {
		  bdp.ssLc.print( fpout ); 
		  fclose( fpout );
	       }
	    }
	    printw("\n> " );
	    getch( );
	    break;

	 default :
	    break;
	 }
      }
      
      nextstar :;

      /*--- plcs and dplcs lightcurves are handled by DPhase class ---*/
      while( plcs[0] )
	 plcs.unhook( plcs[0] );
      while( dplcs[0] )
	 dplcs.unhook( dplcs[0] );

      if( single ) {
	 endwin( );
	 if( fpid ) 
	    fclose( fpid );
	 return 0;
      }
   }
   

   if( fpid ) 
      fclose( fpid );
   endwin( );
}
